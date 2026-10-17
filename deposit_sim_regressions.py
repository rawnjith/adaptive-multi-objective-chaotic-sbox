import inspect
import json
import math
import sys
from pathlib import Path
from argparse import Namespace
from unittest.mock import patch

import numpy as np

TEST_DIR = Path(__file__).resolve().parent
if (TEST_DIR / "sim.py").exists():
    sys.path.insert(0, str(TEST_DIR))
else:
    ROOT = TEST_DIR.parent
    sys.path.insert(0, str(ROOT / "tools"))
import sim  # noqa: E402


def test_lyapunov_uncoupled_mu4_matches_log_two():
    assert math.isclose(sim.lyapunov_2d_logistic(.314159,.271828,4.,0.,n=20_000,transient=2_000), math.log(2), abs_tol=.03)


def test_main_accepts_seed_fpa_and_figure_configuration():
    assert {"seed","fpa_enabled","generate_figures"}.issubset(inspect.signature(sim.main).parameters)


def test_runtime_is_not_an_objective():
    a={"le":1.,"nl_min":100,"du":8,"runtime_ms":99.}; b={"le":1.,"nl_min":100,"du":8,"runtime_ms":.1}
    assert not sim.dominates(a,b) and not sim.dominates(b,a)
    assert "runtime_ms" not in sim.objective_keys()


def test_map_step_is_simultaneous_modulo_one():
    assert np.allclose(sim.step_2d_logistic(.2,.7,4.,.5),((4*.2*.8+.5*.7)%1,(4*.7*.3+.5*.2)%1))


def test_seed_no_fpa_cli_records_config():
    with patch.object(sim.argparse.ArgumentParser,"parse_args",return_value=Namespace(seed=17,no_fpa=True,pop_size=4,generations=1,fpa_interval=1,output="tiny.json",figures=False)):
        with patch.object(sim,"main") as run: sim.cli_main()
    assert run.call_args.kwargs["seed"]==17 and run.call_args.kwargs["fpa_enabled"] is False


def test_main_records_seed_and_fpa_status(tmp_path):
    path=tmp_path/"tiny.json"; sim.main(seed=29,fpa_enabled=False,pop_size=4,n_gen=1,output_path=str(path),generate_figures=False)
    d=json.loads(path.read_text()); assert d["provenance"]["seed"]==29 and d["provenance"]["fpa_enabled"] is False
    assert d["fpa_interval"]>d["n_gen"]


def test_knee_selection_uses_normalized_distance_to_ideal():
    c=[{"le":.8,"nl_min":100,"du":10},{"le":.7,"nl_min":104,"du":10},{"le":.6,"nl_min":98,"du":12}]
    assert sim.select_knee(c)==1


def test_tournament_prefers_rank_then_crowding():
    assert sim.tournament_winner(0,1,[1,0,0],[100,.1,5])==1
    assert sim.tournament_winner(1,2,[1,0,0],[100,.1,5])==2


def test_fisher_yates_steps_simultaneously():
    states=[]; original=sim.step_2d_logistic
    def record(*args): states.append(args[:2]); return original(*args)
    with patch.object(sim,"step_2d_logistic",side_effect=record): sim.map_to_sbox(.3,.7,3.8,.2,transient=2,target_len=4)
    assert len(states)==5


def test_lap_returns_probability_bias():
    assert math.isclose(sim.lap(np.arange(256)),.5,abs_tol=1e-12)


def test_orbit_returns_n_post_transient_states():
    xs,ys=sim.iterate_2d_logistic(.2,.7,4,.5,n=3,transient=0); x,y=.2,.7; expected=[]
    for _ in range(3): x,y=sim.step_2d_logistic(x,y,4,.5); expected.append((x,y))
    assert np.allclose(xs,[z[0] for z in expected]) and np.allclose(ys,[z[1] for z in expected])


def test_periodic_gene_wrapping():
    assert np.allclose(sim.wrap_parameters([1.2,-.2,4.2,.6]),[.2,.8,3.77,.1])


def test_initial_state_transform_is_in_unit_square():
    x,y=sim.initial_state(0,.25); assert 0<=x<1 and 0<=y<1
    assert math.isclose(x,0,abs_tol=1e-12) and math.isclose(y,.5,abs_tol=1e-12)


def test_aes_sbox_generator_matches_known_prefix_and_is_permutation():
    aes=sim.aes_sbox(); assert aes[:4].tolist()==[0x63,0x7c,0x77,0x7b]; assert len(np.unique(aes))==256


def test_aes_sbox_metrics_are_standard_reference_values():
    aes=sim.aes_sbox(); assert sim.nonlinearity(aes)[0]==112; assert sim.differential_uniformity(aes)==4
    assert math.isclose(sim.lap(aes),.0625,abs_tol=1e-12)
