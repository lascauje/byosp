"""Unit tests for Convex-Concave Backtracking for Inertial \
Bregman Proximal Gradient Algorithms in Non-Convex Optimization."""


from bpg import bpg
import numpy as np
from pytest import approx


def test_f() -> None:
    """Sanity check."""
    assert bpg.f(0.1, np.array([[1], [1]])) == approx(0.2)


def test_g() -> None:
    """Sanity check."""
    assert bpg.g([np.array([[1, 3], [2, 4]]), np.array([[1, 3], [2, 4]])],
                 np.array([[1], [2]]),
                 [0.60276, 0.64589]) == approx(347.838348377425)


def test_d_g() -> None:
    """Sanity check."""
    assert sum(bpg.d_g([np.array([[1, 3], [2, 4]]),
                        np.array([[1, 3], [2, 4]])],
                       np.array([[1], [2]]),
                       [0.60276, 0.64589])) == approx(896.77295)


def test_h() -> None:
    """Sanity check."""
    assert bpg.h(np.array([[1], [2]])) == approx(8.75)


def test_d_h() -> None:
    """Sanity check."""
    assert sum(bpg.d_h(np.array([[1], [2]]))) == approx(18.0)


def test_obj_func() -> None:
    """Sanity check."""
    assert bpg.obj_func([np.array([[1, 3], [2, 4]]),
                         np.array([[1, 3], [2, 4]])],
                        np.array([[1], [2]]),
                        [0.60276, 0.64589],
                        0.1) == approx(348.138348377425)


def test_breg_dist() -> None:
    """Sanity check."""
    assert bpg.breg_dist(np.array([[3], [4]]),
                         np.array([[1], [2]])) == approx(124.0)


def test_find_y_test() -> None:
    """Sanity check."""
    assert sum(bpg.find_y(0.9, np.array([[3], [4]]),
                          np.array([[10], [20]]))) == approx(50.7)


def test_find_gamma() -> None:
    """Sanity check."""
    y, gamma = bpg.find_gamma(0.15, 0.00001, np.array([[3], [4]]),
                              np.array([[10], [20]]), 0.001, 10, 0.1)
    assert sum(y) == approx(35.846291405156705)
    assert gamma == approx(0.2541865828329001)


def test_soft_thresh() -> None:
    """Sanity check."""
    assert sum(bpg.soft_thresh(0.01, np.array([[-1.92486], [0.15026]]))) \
        == approx(-1.7746)


def test_t_star() -> None:
    """Sanity check."""
    assert bpg.t_star(np.array([[-1.92486], [0.15026]])) \
        == approx(0.5088523346606675)


def test_find_x() -> None:
    """Sanity check."""
    assert sum(bpg.find_x(0.1,
                          [np.array([[1, 3], [2, 4]]),
                           np.array([[1, 3], [2, 4]])],
                          [0.60276, 0.64589],
                          np.array([[-1], [0]]),
                          10)) == approx(-0.9054999300267382)


def test_lyapunov() -> None:
    """Sanity check."""
    assert bpg.lyapunov(0.1, [np.array([[1, 3], [2, 4]]),
                              np.array([[1, 3], [2, 4]])],
                        [0.60276, 0.64589],
                        np.array([[3], [4]]),
                        np.array([[1], [2]]),
                        10) == approx(238.8138348377425)


def test_find_lower() -> None:
    """Sanity check."""
    lower, y, gamma = bpg.find_lower(0.15,
                                     0.00001,
                                     [np.array([[0.99531367, 0.58327082],
                                                [0.9877355, 0.66336516]]),
                                      np.array([[0.0981552, 0.52265855],
                                                [0.61172762, 0.0371408]])],
                                     [0.6027633760716439, 0.6458941130666561],
                                     np.array([[0.4719454622448402],
                                               [0.4220473240632225]]),
                                     np.array([[0.4693035792040524],
                                               [0.42080845789332616]]),
                                     0.001,
                                     10)
    assert lower == approx(0.512)
    assert sum(y) == approx(0.8887589035161779)
    assert gamma == approx(0.34867844010000015)


def test_find_upper() -> None:
    """Sanity check."""
    upper, x = bpg.find_upper(0.1,
                              [np.array([[0.99531367, 0.58327082],
                                         [0.9877355, 0.66336516]]),
                               np.array([[0.0981552, 0.52265855],
                                         [0.61172762, 0.0371408]])],
                              [0.6027633760716439, 0.6458941130666561],
                              np.array([[0.4588713468854095],
                                        [0.4191807008769411]]),
                              32.768,
                              10)

    assert upper == approx(80.0)
    assert sum(x) == approx(0.8779832567675403)


def test_find_upper_official() -> None:
    """Sanity check."""
    upper, x = bpg.find_upper_official(0.1,
                                       [np.array([[0.99531367, 0.58327082],
                                                  [0.9877355, 0.66336516]]),
                                        np.array([[0.0981552, 0.52265855],
                                                  [0.61172762, 0.0371408]])],
                                       [0.6027633760716439,
                                        0.6458941130666561],
                                       np.array([[0.4588713468854095],
                                                 [0.4191807008769411]]),
                                       32.768,
                                       10)

    assert upper == approx(160.0)
    assert sum(x) == approx(0.8780176535863605)


def test_lcg() -> None:
    """Sanity check."""
    assert sum(bpg.lcg(10, 0)) == approx(3.977038038428873)


def test_iter_coca_bpg_10() -> None:
    """Sanity check.

    Lyapunov values must be decreasing.
    """
    dim = 10
    nb_time = 100
    res = \
        bpg.iter_coca_bpg(lam=0.1,
                          delta=0.15,
                          epsilon=0.00001,
                          a=bpg.init_a(nb_time, dim),
                          b=bpg.init_b(nb_time),
                          prev_x=np.ones(dim),
                          x=np.ones(dim),
                          lower=0.001,
                          upper=10,
                          nb_iter=1000)

    lyapunov_vals = res[-1]
    assert all(lyapunov_vals[i]
               >= lyapunov_vals[i + 1] for i in range(len(lyapunov_vals) - 1)) # noqa
    assert lyapunov_vals[-1] == approx(0.0009734142817317467)


def test_iter_coca_bpg_128() -> None:
    """Sanity check.

    Lyapunov values must be decreasing.
    """
    dim = 128
    nb_time = 100
    res = \
        bpg.iter_coca_bpg(lam=0.1,
                          delta=0.15,
                          epsilon=0.00001,
                          a=bpg.init_a(nb_time, dim),
                          b=bpg.init_b(nb_time),
                          prev_x=np.ones(dim),
                          x=np.ones(dim),
                          lower=0.001,
                          upper=10,
                          nb_iter=1000)

    lyapunov_vals = res[-1]
    assert all(lyapunov_vals[i]
               >= lyapunov_vals[i + 1] for i in range(len(lyapunov_vals) - 1)) # noqa
    assert lyapunov_vals[-1] == approx(8.698883636567947e-05)


def test_official() -> None:
    """Sanity check for the official solution."""
    lyapunov_vals = bpg.main_official()
    assert all(lyapunov_vals[i]
               >= lyapunov_vals[i + 1] for i in range(len(lyapunov_vals) - 1)) # noqa
    assert lyapunov_vals[-1] == approx(0.0010232284321974865)
