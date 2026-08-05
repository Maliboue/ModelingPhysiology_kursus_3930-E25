
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from scipy.integrate import solve_ivp


# Noble (1962) Purkinje-fibre model parameters.
# Units:
#   time: ms
#   voltage: mV
#   current density: microA / cm^2
#   conductance density: mmho / cm^2
#   capacitance: microF / cm^2
params = namedtuple(
    "params",
    [
        "C",
        "g_Na_max",
        "g_Na_background",
        "g_K1_scale",
        "g_K2_max",
        "E_Na",
        "E_K",
        "g_An",
        "E_An",
        "I_app",
    ],
)


def _x_over_expm1(x):
    """Stable evaluation of x / (exp(x) - 1), including x = 0."""
    x = np.asarray(x)
    small = np.abs(x) < 1e-7
    result = np.empty_like(x, dtype=float)
    result[small] = 1.0 - x[small] / 2.0 + x[small] ** 2 / 12.0
    result[~small] = x[~small] / np.expm1(x[~small])
    return float(result) if result.ndim == 0 else result


def noble_1962_model(t, state, model_parameters):
    """
    Noble (1962) Purkinje-fibre membrane model, equations (1)-(20).

    Parameters
    ----------
    t : float
        Time in ms.

    state : tuple with 4 elements
        Current values of:
            E_m : membrane potential in mV
            m   : sodium activation gate
            h   : sodium inactivation gate
            n   : delayed potassium activation gate

    model_parameters : collections.namedtuple
        Model parameters stored as attributes.

    Returns
    -------
    tuple with 4 elements
        dE_m/dt, dm/dt, dh/dt, dn/dt.
    """

    C = model_parameters.C

    g_Na_max = model_parameters.g_Na_max
    g_Na_background = model_parameters.g_Na_background
    g_K1_scale = model_parameters.g_K1_scale
    g_K2_max = model_parameters.g_K2_max

    E_Na = model_parameters.E_Na
    E_K = model_parameters.E_K

    g_An = model_parameters.g_An
    E_An = model_parameters.E_An

    I_app = model_parameters.I_app

    E_m = state[0]
    m = state[1]
    h = state[2]
    n = state[3]

    # Potassium rate constants: equations (8)-(9)
    alpha_n = 0.001 * _x_over_expm1((-E_m - 50.0) / 10.0)
    beta_n = 0.002 * np.exp((-E_m - 90.0) / 80.0)

    # Sodium h-gate rate constants: equations (16)-(17)
    alpha_h = 0.17 * np.exp((-E_m - 90.0) / 20.0)
    beta_h = 1.0 / (np.exp((-E_m - 42.0) / 10.0) + 1.0)

    # Sodium m-gate rate constants: equations (18)-(19)
    alpha_m = 1.5 * _x_over_expm1((-E_m - 48.0) / 15.0)
    beta_m = 0.6 * _x_over_expm1((E_m + 8.0) / 5.0)

    # Instantaneous and delayed potassium conductances: equations (5)-(6)
    g_K1 = g_K1_scale * (
        np.exp((-E_m - 90.0) / 50.0)
        + 0.0125 * np.exp((E_m + 90.0) / 60.0)
    )
    g_K2 = g_K2_max * n**4

    # Ionic currents: equations (3), (10), and (20)
    I_Na = (g_Na_max * m**3 * h + g_Na_background) * (E_m - E_Na)
    I_K = (g_K1 + g_K2) * (E_m - E_K)
    I_An = g_An * (E_m - E_An)

    # Membrane equation (4), using outward-positive current convention.
    dE_m_dt = (I_app - I_Na - I_K - I_An) / C

    # Gating equations (7), (13), and (14)
    dm_dt = alpha_m * (1.0 - m) - beta_m * m
    dh_dt = alpha_h * (1.0 - h) - beta_h * h
    dn_dt = alpha_n * (1.0 - n) - beta_n * n

    return (dE_m_dt, dm_dt, dh_dt, dn_dt)


def steady_state_gates(E_m):
    """Return m_inf, h_inf, and n_inf at a specified membrane potential."""
    alpha_n = 0.001 * _x_over_expm1((-E_m - 50.0) / 10.0)
    beta_n = 0.002 * np.exp((-E_m - 90.0) / 80.0)

    alpha_h = 0.17 * np.exp((-E_m - 90.0) / 20.0)
    beta_h = 1.0 / (np.exp((-E_m - 42.0) / 10.0) + 1.0)

    alpha_m = 1.5 * _x_over_expm1((-E_m - 48.0) / 15.0)
    beta_m = 0.6 * _x_over_expm1((E_m + 8.0) / 5.0)

    m_inf = alpha_m / (alpha_m + beta_m)
    h_inf = alpha_h / (alpha_h + beta_h)
    n_inf = alpha_n / (alpha_n + beta_n)

    return m_inf, h_inf, n_inf


# Parameters used for the spontaneous rhythm in Figure 6A.
# Figure 6A has no anion current, so g_An = 0.
model_parameters = params(
    C=12.0,
    g_Na_max=400.0,
    g_Na_background=0.14,
    g_K1_scale=1.2,
    g_K2_max=1.2,
    E_Na=40.0,
    E_K=-100.0,
    g_An=0.0,
    E_An=-60.0,   # irrelevant while g_An = 0
    I_app=0.0,
)


# Start near the middle of the pacemaker potential, as described in the paper.
E_m0 = -80.0
m0, h0, n_inf0 = steady_state_gates(E_m0)

# Noble states that the initial n value was guessed and loses influence after
# one cycle. A value of 0.30 gives rapid convergence to the stable rhythm.
initial_state = (E_m0, m0, h0, 0.30)


# First integrate long enough to remove dependence on the guessed initial n.
burn_in = solve_ivp(
    fun=noble_1962_model,
    t_span=(0.0, 5000.0),
    y0=initial_state,
    args=(model_parameters,),
    method="LSODA",
    rtol=1e-9,
    atol=1e-11,
    max_step=0.1,
)

if not burn_in.success:
    raise RuntimeError(burn_in.message)


# Reproduce approximately two cycles, matching the 0-2 s extent of Fig. 6A.
plot_duration = 2000.0
t_eval = np.arange(0.0, plot_duration + 0.1, 0.1)

solution = solve_ivp(
    fun=noble_1962_model,
    t_span=(0.0, plot_duration),
    y0=burn_in.y[:, -1],
    t_eval=t_eval,
    args=(model_parameters,),
    method="LSODA",
    rtol=1e-9,
    atol=1e-11,
    max_step=0.1,
)

if not solution.success:
    raise RuntimeError(solution.message)


fig, ax = plt.subplots(figsize=(8.0, 4.5))
ax.plot(solution.t / 1000.0, solution.y[0], linewidth=1.5)
ax.axhline(0.0, linewidth=0.7)
ax.set_xlim(0.0, 2.0)
ax.set_ylim(-100.0, 50.0)
ax.set_xlabel("Time (s)")
ax.set_ylabel(r"$E_m$ (mV)")
ax.set_title("Noble (1962) Purkinje-fibre model - reproduction of Fig. 6A")
fig.tight_layout()
plt.show()
