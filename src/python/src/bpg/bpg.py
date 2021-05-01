"""Implementation of Convex-Concave Backtracking for Inertial \
Bregman Proximal Gradient Algorithms in Non-Convex Optimization."""


from timeit import default_timer as timer
from typing import Callable, Tuple

import numpy as np


def f(lam: float, x: np.ndarray) -> float:
    """See p.20 (6.5)."""
    return lam * np.linalg.norm(x, ord=1)


def g(a: list, x: np.ndarray, b: list) -> float:
    """See p.19 (6.5)."""
    def g_helper(a: np.ndarray, x: np.ndarray, b: float) -> float:
        return 0.25 * ((np.dot(x.T, np.dot(a, x)) - b) ** 2)
    return sum([g_helper(a[i], x, b[i]) for i in range(len(a))])


def d_g(a: list, x: np.ndarray, b: list) -> np.ndarray:
    """See p.19 (6.5)."""
    def d_g_helper(a: np.ndarray, x: np.ndarray, b: float) -> float:
        # It should be "b**2, but the official implementation uses "b"
        return (np.dot(x.T, np.dot(a, x)) - b) * np.dot(a, x)
    return sum([d_g_helper(a[i], x, b[i]) for i in range(len(a))])


def h(x: np.ndarray) -> float:
    """See p.20 (6.6)."""
    return 0.25 * (np.linalg.norm(x) ** 4) + 0.5 * (np.linalg.norm(x) ** 2)


def d_h(x: np.ndarray) -> np.ndarray:
    """See p.20 (6.6)."""
    return np.dot(x.T, x) * x + x


def obj_func(a: list, x: np.ndarray, b: list, lam: float) -> float:
    """See p.19 (6.5)."""
    return f(lam, x) + g(a, x, b)


def breg_dist(x: np.ndarray, y: np.ndarray) -> float:
    """See p.4 definition 2.1."""
    dist = h(x) - h(y) - np.dot((x - y).T, d_h(y))
    return dist if 1e-15 < dist else 0


def find_y(gamma: float, prev_x: np.ndarray, x: np.ndarray) -> np.ndarray:
    """See p.7 (3.4)."""
    return x + gamma * (x - prev_x)


def find_gamma(delta: float, epsilon: float,
               prev_x: np.ndarray, x: np.ndarray,
               lower: float, upper: float,
               prev_tau: float) -> Tuple[np.ndarray, float]:
    """See p.7 (3.5)."""
    gamma = 1
    y = find_y(gamma, prev_x, x)
    while(((delta - epsilon) * breg_dist(prev_x, x) / (1 + (lower * prev_tau)))
          < breg_dist(x, y)): # noqa
        gamma = gamma * 0.9
        y = find_y(gamma, prev_x, x)
    return y, gamma


def find_lower(delta: float, epsilon: float,
               a: list, b: list,
               prev_x: np.ndarray, x: np.ndarray,
               lower: float, upper: float) -> Tuple[float, np.ndarray, float]:
    """See p.7 (3.5)."""
    while True:
        y, gamma = \
            find_gamma(delta, epsilon, prev_x, x, lower, upper, 1 / upper)
        if not 1e-7 \
           < g(a, y, b) + np.dot(d_g(a, y, b).T, x - y) \
           - lower * breg_dist(x, y) - g(a, x, b):
            break
        else:
            # Scaling "lower L" parameter with 2
            lower = lower * 2
    return lower, y, gamma


def soft_thresh(theta: float, y: np.ndarray) -> np.ndarray:
    """See p.21 (6.10)."""
    return np.maximum((np.abs(y) - theta), 0) * np.sign(y)


def t_star(coeff: np.ndarray) -> float:
    """See p.21 (6.10) next."""
    return np.roots([np.linalg.norm(coeff) ** 2, 0, 1, -1])[-1].real


def find_x(lam: float,
           a: list, b: list,
           y: np.ndarray,
           upper: float) -> np.ndarray:
    """See both p.7 (3.7) and p.21 (6.10) next."""
    tau = 1 / upper
    s = soft_thresh(lam * tau, d_h(y) - tau * d_g(a, y, b))
    return t_star(s) * s


def find_upper(lam: float,
               a: list, b: list, y: np.ndarray,
               lower: float, upper: float) -> Tuple[float, np.ndarray]:
    """See p.7 (3.5)."""
    while True:
        x = find_x(lam, a, b, y, upper)
        if not g(a, y, b) + np.dot(d_g(a, y, b).T, (x - y)) \
           + upper * breg_dist(x, y) - g(a, x, b) < -1e-7:
            break
        else:
            # Scaling "upper L" parameter with 2
            upper = upper * 2
    return upper, x


def find_upper_official(lam: float,
                        a: list, b: list,
                        y: np.ndarray,
                        lower: float, upper: float) -> Tuple[float,
                                                             np.ndarray]:
    """See p.7 (3.5).

    This is a wrapper to the raw find_upper function.
    In the raw find_upper function the paper has been implemented
    slightly differently, but in order to check the results
    this function (find_upper_official) is used instead.
    """
    x = find_x(lam, a, b, y, upper)
    if g(a, y, b) + np.dot(d_g(a, y, b).T, (x - y)) \
       + upper * breg_dist(x, y) - g(a, x, b) < -1e-7:
        upper, _ = find_upper(lam, a, b, y, lower, upper)
        next_upper = upper * 2
        next_x = find_x(lam, a, b, y, next_upper)
        return next_upper, next_x
    else:
        return upper, x


def lyapunov(lam: float,
             a: list, b: list,
             prev_x: np.ndarray, x: np.ndarray,
             upper: float) -> float:
    """See p.11 (5.3) (simpler)."""
    prev_tau = 1 / upper
    return prev_tau * obj_func(a, x, b, lam) + breg_dist(x, prev_x)


def coca_bpg(lam: float, delta: float, epsilon: float,
             a: list, b: list,
             prev_x: np.ndarray, x: np.ndarray,
             prev_lower: float, prev_upper: float) -> Tuple[
                 float,
                 np.ndarray, np.ndarray, np.ndarray,
                 float, float, float]:
    """See p.7 Convex-Concave Inertial BPG algorithm."""
    lower, y, gamma = \
        find_lower(delta, epsilon, a, b, prev_x, x, prev_lower, prev_upper)
    upper, next_x = find_upper_official(lam, a, b, y, lower, prev_upper)
    res = lyapunov(lam, a, b, x, next_x, upper)
    return gamma, x, next_x, y, lower, upper, res


def iter_coca_bpg(lam: float, delta: float, epsilon: float,
                  a: list, b: list, prev_x: np.ndarray, x: np.ndarray,
                  lower: float, upper: float,
                  nb_iter: int) -> Tuple[float, float, float,
                                         list, list,
                                         np.ndarray, np.ndarray,
                                         float, float, list]:
    """Wrapper to iterate on coca-bpg function."""
    res_all = []
    for _i in range(nb_iter):
        gamma, prev_x, x, y, lower, upper, res = \
            coca_bpg(lam, delta, epsilon, a, b, prev_x, x, lower, upper)
        res_all = res_all + [res]
    return lam, delta, epsilon, a, b, prev_x, x, lower, upper, res_all


def lcg(nb_time: int, seed: int) -> list:
    """A simple linear congruential generator inspired from the glibc.

    Values are in range ]0,1[.
    See glibc/stdlib/random_r.c:__random_r:341
    """

    def lcg_glibc(x: int) -> Callable[[], int]:
        def rand() -> int:
            nonlocal x
            x = (1103515245 * x + 12345) & 0x7fffffff
            return x
        return rand
    values = [seed]
    for _i in range(nb_time):
        values.append(lcg_glibc(values[-1])())
    return [v / 0x7fffffff for v in values][1:]


def init_a(nb_time: int, row: int) -> list:
    """Init of a (list of matrices) using a linear congruential generator.

    See p.19 (6.3).
    """
    def single_matrix(seed: int) -> np.ndarray:
        m = np.array(lcg(row, seed)).reshape(row, 1)
        return np.dot(m, m.T)
    return [*map(single_matrix, range(nb_time))]


def init_b(nb_time: int) -> list:
    """Init of b (list of float) using a linear congruential generator.

    See p.19 (6.3).
    """
    return lcg(nb_time, 0)


def bench(dim: int, nb_time: int = 100, nb_iter: int = 1000) -> None:
    """Benchmark the BPG algorithm."""
    a = init_a(nb_time, dim)
    b = init_b(nb_time)
    lam = 1e-1
    delta = 0.15
    epsilon = 0.00001
    upper = 10
    lower = 1e-4 * upper
    x = np.ones(dim)
    prev_x = x
    start = timer()
    lyapunov_vals = \
        iter_coca_bpg(lam, delta, epsilon, a, b, prev_x, x,
                      lower, upper, nb_iter)[-1]
    end = timer()
    print(str((end - start) * 1000) + " msecs")
    return lyapunov_vals[-1]


def main_official() -> list:
    """The official usage of the BPG algorithm as described in the paper."""
    np.random.seed(0)
    nb_iter = 1000
    dim = 10
    temp_Alist = []
    temp_blist = []

    for _i in range(100):
        temp_A = np.random.rand(dim, 1)
        A = temp_A * temp_A.T
        temp_Alist = temp_Alist + [A]
        temp_b = np.random.rand(1)[0]
        temp_blist = temp_blist + [temp_b]

    a = temp_Alist
    b = temp_blist

    lam = 1e-1
    delta = 0.15
    epsilon = 0.00001
    upper = 10
    lower = 1e-4 * upper
    x = np.ones(dim)
    prev_x = x
    start = timer()
    lyapunov_vals = \
        iter_coca_bpg(lam, delta, epsilon, a, b, prev_x, x,
                      lower, upper, nb_iter)[-1]
    end = timer()
    print(str((end - start) * 1000) + " msecs")
    print("lyapunov_vals: ")
    print(lyapunov_vals[-1])
    print(all(lyapunov_vals[i]
              >= lyapunov_vals[i + 1] for i in range(len(lyapunov_vals) - 1))) # noqa

    return lyapunov_vals


def main() -> None:
    """I don't do a whole lot ... yet."""
    print("Hello BPG!")
