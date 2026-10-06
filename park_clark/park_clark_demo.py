"""Clarke/Park-Transformation fuer Umrichter und MMC-HVDC - Demo zum Lernskript.

Ausfuehren:  python3 park_clark/park_clark_demo.py
Benoetigt numpy, optional matplotlib (fuer Plots in park_clark/plots/).

Konvention: amplitudeninvariant (Faktor 2/3), cos-basiert, d-Achse auf Phase a
bei theta = 0, q-Achse 90 Grad voraus.  x_dq = x_alphabeta * exp(-j*theta).
"""

from pathlib import Path

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:  # Plots sind optional
    plt = None

F_GRID = 50.0
OMEGA = 2 * np.pi * F_GRID
PLOT_DIR = Path(__file__).parent / "plots"


# --------------------------------------------------------------------------
# Transformationen
# --------------------------------------------------------------------------
def clarke(a, b, c):
    """abc -> (alpha, beta, zero), amplitudeninvariant."""
    alpha = 2 / 3 * (a - 0.5 * b - 0.5 * c)
    beta = (b - c) / np.sqrt(3)
    zero = (a + b + c) / 3
    return alpha, beta, zero


def inv_clarke(alpha, beta, zero=0.0):
    a = alpha + zero
    b = -0.5 * alpha + np.sqrt(3) / 2 * beta + zero
    c = -0.5 * alpha - np.sqrt(3) / 2 * beta + zero
    return a, b, c


def park(alpha, beta, theta):
    """(alpha, beta) -> (d, q):  d + jq = (alpha + j beta) e^{-j theta}."""
    d = alpha * np.cos(theta) + beta * np.sin(theta)
    q = -alpha * np.sin(theta) + beta * np.cos(theta)
    return d, q


def inv_park(d, q, theta):
    alpha = d * np.cos(theta) - q * np.sin(theta)
    beta = d * np.sin(theta) + q * np.cos(theta)
    return alpha, beta


def abc_to_dq(a, b, c, theta):
    alpha, beta, zero = clarke(a, b, c)
    d, q = park(alpha, beta, theta)
    return d, q, zero


def three_phase(amp, phi, t, omega=OMEGA, seq=+1):
    """Symmetrisches Drehstromsystem. seq=+1 Mitsystem (a-b-c), -1 Gegensystem."""
    return (
        amp * np.cos(omega * t + phi),
        amp * np.cos(omega * t + phi - seq * 2 * np.pi / 3),
        amp * np.cos(omega * t + phi + seq * 2 * np.pi / 3),
    )


def save(fig, name):
    PLOT_DIR.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / name, dpi=120)
    plt.close(fig)
    print(f"   Plot gespeichert: {PLOT_DIR / name}")


# --------------------------------------------------------------------------
# 1) Grundlagen: abc -> alpha/beta -> dq
# --------------------------------------------------------------------------
def demo_basics():
    print("\n1) Clarke + Park auf ein symmetrisches System")
    t = np.linspace(0, 0.04, 2000)
    theta = OMEGA * t
    va, vb, vc = three_phase(100.0, np.deg2rad(30), t)
    alpha, beta, zero = clarke(va, vb, vc)
    vd, vq, v0 = abc_to_dq(va, vb, vc, theta)

    print(f"   v_d = {vd.mean():7.3f} (Schwankung {np.ptp(vd):.1e})")
    print(f"   v_q = {vq.mean():7.3f} (Schwankung {np.ptp(vq):.1e})")
    print(f"   v_0 = {v0.mean():7.3f}")
    print("   -> Uebung 1: 100*cos(30deg)=86.603, 100*sin(30deg)=50.000")
    assert np.allclose(vd, 100 * np.cos(np.pi / 6))
    assert np.allclose(vq, 100 * np.sin(np.pi / 6))

    # Ruecktransformation muss wieder abc liefern
    a2, b2, c2 = inv_clarke(*inv_park(vd, vq, theta), v0)
    assert np.allclose([a2, b2, c2], [va, vb, vc])
    print("   Ruecktransformation dq -> abc: OK")

    if plt:
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
        ax[0].plot(t * 1e3, va, t * 1e3, vb, t * 1e3, vc)
        ax[0].set(title="abc", xlabel="t [ms]")
        ax[0].legend(["a", "b", "c"])
        ax[1].plot(t * 1e3, alpha, t * 1e3, beta)
        ax[1].set(title="alpha/beta (Clarke)", xlabel="t [ms]")
        ax[1].legend(["alpha", "beta"])
        ax[2].plot(t * 1e3, vd, t * 1e3, vq)
        ax[2].set(title="dq (Park) -> Gleichgroessen", xlabel="t [ms]")
        ax[2].legend(["d", "q"])
        for a in ax:
            a.grid(alpha=0.3)
        save(fig, "1_abc_alphabeta_dq.png")


# --------------------------------------------------------------------------
# 2) Leistung: abc vs. dq
# --------------------------------------------------------------------------
def demo_power():
    print("\n2) Leistung in dq")
    t = np.linspace(0, 0.02, 1000)
    theta = OMEGA * t
    V, I, phi_i = 100.0, 10.0, np.deg2rad(-25)  # Strom eilt 25 Grad nach
    va, vb, vc = three_phase(V, 0.0, t)
    ia, ib, ic = three_phase(I, phi_i, t)

    p_abc = va * ia + vb * ib + vc * ic
    # Blindleistung aus abc ueber verkettete Spannungen
    q_abc = ((vb - vc) * ia + (vc - va) * ib + (va - vb) * ic) / np.sqrt(3)

    vd, vq, _ = abc_to_dq(va, vb, vc, theta)
    id_, iq, _ = abc_to_dq(ia, ib, ic, theta)
    p_dq = 1.5 * (vd * id_ + vq * iq)
    q_dq = 1.5 * (vq * id_ - vd * iq)

    p_ref = 1.5 * V * I * np.cos(-phi_i)
    q_ref = 1.5 * V * I * np.sin(-phi_i)
    print(f"   P: abc={p_abc.mean():8.2f}  dq={p_dq.mean():8.2f}  Formel={p_ref:8.2f}")
    print(f"   Q: abc={q_abc.mean():8.2f}  dq={q_dq.mean():8.2f}  Formel={q_ref:8.2f}")
    print(f"   i_d={id_.mean():.3f}, i_q={iq.mean():.3f}  (v_q=0 -> P~i_d, Q~-i_q)")
    assert np.allclose(p_abc, p_dq) and np.allclose(q_abc, q_dq)


# --------------------------------------------------------------------------
# 3) SRF-PLL
# --------------------------------------------------------------------------
def demo_pll():
    print("\n3) SRF-PLL: Frequenz 50.5 Hz und Phasensprung 30 Grad bei t=0.1 s")
    dt, t_end = 1e-5, 0.25
    t = np.arange(0, t_end, dt)
    f_true = 50.5
    jump = np.where(t >= 0.1, np.deg2rad(30), 0.0)
    va, vb, vc = three_phase(1.0, jump, t, omega=2 * np.pi * f_true)
    alpha, beta, _ = clarke(va, vb, vc)

    # PI-Auslegung fuer v_d ~ 1 p.u.:  Kp = 2*zeta*wn, Ki = wn^2
    wn, zeta = 2 * np.pi * 20, 0.707
    kp, ki = 2 * zeta * wn, wn**2
    theta_hat, integ = 0.0, 0.0
    f_hat = np.zeros_like(t)
    vq_log = np.zeros_like(t)
    for k in range(len(t)):
        _, vq = park(alpha[k], beta[k], theta_hat)
        integ += ki * vq * dt
        w_hat = OMEGA + kp * vq + integ
        theta_hat = (theta_hat + w_hat * dt) % (2 * np.pi)
        f_hat[k] = w_hat / (2 * np.pi)
        vq_log[k] = vq
    print(f"   geschaetzte Frequenz am Ende: {f_hat[-1]:.3f} Hz (wahr {f_true} Hz)")
    print(f"   v_q am Ende: {vq_log[-1]:.1e} -> d-Achse liegt auf dem Spannungszeiger")

    if plt:
        fig, ax = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
        ax[0].plot(t, f_hat)
        ax[0].axhline(f_true, ls="--", c="k")
        ax[0].set(ylabel="f_hat [Hz]", title="SRF-PLL")
        ax[1].plot(t, vq_log)
        ax[1].set(ylabel="v_q [p.u.]", xlabel="t [s]")
        for a in ax:
            a.grid(alpha=0.3)
        save(fig, "3_pll.png")


# --------------------------------------------------------------------------
# 4) dq-Stromregelung mit und ohne Entkopplung
# --------------------------------------------------------------------------
def simulate_current_control(decouple, L=5e-3, R=0.05, Vg=325.0, dt=1e-5, t_end=0.06):
    """Gemitteltes Modell in dq: L di/dt = v_c - v_g - R i - j w L i."""
    alpha_c = 2 * np.pi * 300  # Bandbreite des Stromreglers
    kp, ki = alpha_c * L, alpha_c * R
    t = np.arange(0, t_end, dt)
    i = 0j
    integ = 0j
    log = np.zeros(len(t), complex)
    for k, tk in enumerate(t):
        i_ref = (100.0 if tk >= 0.01 else 0.0) + 1j * (-50.0 if tk >= 0.035 else 0.0)
        err = i_ref - i
        integ += ki * err * dt
        v_c = Vg + kp * err + integ  # Netzspannungsvorsteuerung (v_gd=Vg, v_gq=0)
        if decouple:
            v_c += 1j * OMEGA * L * i  # d: -wL i_q, q: +wL i_d
        di = (v_c - Vg - R * i - 1j * OMEGA * L * i) / L
        i += di * dt
        log[k] = i
    return t, log


def demo_current_control():
    print("\n4) dq-Stromregler: Sprung i_d* = 100 A (t=10 ms), i_q* = -50 A (t=35 ms)")
    t, with_dec = simulate_current_control(True)
    _, without = simulate_current_control(False)
    win = (t > 0.01) & (t < 0.035)
    print(f"   max |i_q| waehrend i_d-Sprung, mit Entkopplung:  {np.abs(with_dec.imag[win]).max():6.2f} A")
    print(f"   max |i_q| waehrend i_d-Sprung, ohne Entkopplung: {np.abs(without.imag[win]).max():6.2f} A")

    if plt:
        fig, ax = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
        for a, sig, name in [(ax[0], np.real, "i_d"), (ax[1], np.imag, "i_q")]:
            a.plot(t * 1e3, sig(with_dec), label="mit Entkopplung")
            a.plot(t * 1e3, sig(without), "--", label="ohne Entkopplung")
            a.set(ylabel=f"{name} [A]")
            a.grid(alpha=0.3)
            a.legend()
        ax[1].set(xlabel="t [ms]")
        ax[0].set(title="Vektorstromregelung in dq")
        save(fig, "4_current_control.png")


# --------------------------------------------------------------------------
# 5) MMC: Kreisstrom im -2w-Rahmen und Nullsystem
# --------------------------------------------------------------------------
def demo_mmc_circulating():
    print("\n5) MMC-Kreisstrom: DC-Anteil + 2. Harmonische im Gegensystem")
    t = np.linspace(0, 0.04, 4000)
    theta = OMEGA * t
    I_dc = 1500.0
    # i_c,k = I_dc/3 + 2. Harmonische, Gegensystem (Folge a-c-b)
    h2 = three_phase(120.0, np.deg2rad(40), t, omega=2 * OMEGA, seq=-1)
    ic_a, ic_b, ic_c = (I_dc / 3 + x for x in h2)

    alpha, beta, zero = clarke(ic_a, ic_b, ic_c)
    d_pos, q_pos = park(alpha, beta, theta)  # normaler +w-Rahmen
    d_m2, q_m2 = park(alpha, beta, -2 * theta)  # CCSC-Rahmen mit -2 theta

    print(f"   Nullsystem i_0 = {zero.mean():.1f} A = I_dc/3 = {I_dc / 3:.1f} A (Energie DC <-> Arme)")
    print(f"   +w-Rahmen:  i_d schwingt (ptp {np.ptp(d_pos):6.1f} A) -> ungeeignet fuer PI")
    print(f"   -2w-Rahmen: i_d={d_m2.mean():6.1f}, i_q={q_m2.mean():6.1f} A, ptp {np.ptp(d_m2):.1e}")
    print("   -> CCSC: PI-Regler treiben diese Gleichgroessen auf 0")

    if plt:
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
        ax[0].plot(t * 1e3, ic_a, t * 1e3, ic_b, t * 1e3, ic_c)
        ax[0].set(title="Kreisstrom i_c,abc", xlabel="t [ms]")
        ax[1].plot(t * 1e3, d_pos, t * 1e3, q_pos)
        ax[1].set(title="Park mit +theta", xlabel="t [ms]")
        ax[2].plot(t * 1e3, d_m2, t * 1e3, q_m2)
        ax[2].set(title="Park mit -2 theta (CCSC)", xlabel="t [ms]")
        ax[2].legend(["d", "q"])
        for a in ax:
            a.grid(alpha=0.3)
        save(fig, "5_mmc_circulating_current.png")


# --------------------------------------------------------------------------
# 6) Unsymmetrisches Netz: Gegensystem -> 2w-Ripple in dq
# --------------------------------------------------------------------------
def demo_unbalanced():
    print("\n6) Unsymmetrische Netzspannung (10 % Gegensystem)")
    t = np.linspace(0, 0.04, 4000)
    theta = OMEGA * t
    pos = three_phase(1.0, 0.0, t)
    neg = three_phase(0.1, 0.0, t, seq=-1)
    va, vb, vc = (p + n for p, n in zip(pos, neg))
    alpha, beta, _ = clarke(va, vb, vc)
    d_p, q_p = park(alpha, beta, theta)
    d_n, q_n = park(alpha, beta, -theta)
    print(f"   +w-Rahmen: v_d Mittel {d_p.mean():.3f}, Ripple 100 Hz ptp {np.ptp(d_p):.3f}")
    print(f"   -w-Rahmen: v_d Mittel {d_n.mean():.3f} (= Gegensystem), Ripple ptp {np.ptp(d_n):.3f}")
    print("   -> Mittelwerte trennen (Notch 100 Hz / DDSRF) = Mit- und Gegensystem")

    if plt:
        fig, ax = plt.subplots(figsize=(8, 3.6))
        ax.plot(t * 1e3, d_p, label="v_d (+w)")
        ax.plot(t * 1e3, q_p, label="v_q (+w)")
        ax.plot(t * 1e3, d_n, "--", label="v_d (-w)")
        ax.set(title="Unsymmetrie: 100-Hz-Ripple in dq", xlabel="t [ms]")
        ax.grid(alpha=0.3)
        ax.legend()
        save(fig, "6_unbalanced.png")


if __name__ == "__main__":
    demo_basics()
    demo_power()
    demo_pll()
    demo_current_control()
    demo_mmc_circulating()
    demo_unbalanced()
    print("\nAlle Pruefungen bestanden.")
