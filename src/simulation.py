# simulation.py

import numpy as np


# ==========================================================
# Utility functions
# ==========================================================

def gaussian(x, mu, sigma, amp):
    return amp * np.exp(-(x - mu)**2 / (2 * sigma**2))

def gaussian2D(x, y, mu_x, mu_y, sigma_x, sigma_y, amp):
    # Calculate the exponent term
    exponent = -((x - mu_x)**2 / (2 * sigma_x**2) + (y - mu_y)**2 / (2 * sigma_y**2))
    # Calculate the full Gaussian function
    g = amp * np.exp(exponent)
    return g


def normalize(x):
    max_val = np.max(np.abs(x))
    return x if max_val == 0 else x / max_val


# ==========================================================
# 1D Spectral Model
# ==========================================================

def generate_ppm_axis(n_points=3000, start=12, stop=6.5):
    return np.linspace(start, stop, n_points)


def build_spectrum(ppm, peaks):
    spectrum = np.zeros_like(ppm)
    for mu, sigma, amp in peaks:
        spectrum += gaussian(ppm, mu, sigma, amp)
    return spectrum


def generate_pure_spectra(ppm, seed=2):
    """
    Generate three distinct but related 1D spectra:
    holo, intermediate, apo
    """
    np.random.seed(seed)

    # Shared scaffold
    shared_peaks = [
        (11.4,0.020,0.6),(10.9,0.018,0.7),
        (9.6,0.020,0.8),(9.1,0.018,0.6),

        (8.8,0.020,0.9),(8.5,0.018,1.0),
        (8.3,0.020,0.9),(8.1,0.018,0.8),

        (7.9,0.017,1.1),(7.8,0.015,0.9),
        (7.7,0.015,1.0),(7.6,0.016,1.2),
        (7.5,0.014,0.9),(7.3,0.016,1.3),
        (7.2,0.015,0.8),(7.1,0.015,1.0),
        (6.9,0.014,1.5),(6.7,0.015,1.0),
        (6.6,0.015,1.2),(6.55,0.014,0.8),
    ]

    # --- HOLO ---
    holo_peaks = shared_peaks + [
        (8.65,0.016,0.8),
        (7.55,0.013,0.9),
        (6.95,0.013,0.7),
    ]
    holo = build_spectrum(ppm, holo_peaks)

    # --- INTERMEDIATE ---
    int_peaks = []
    for mu, sig, amp in shared_peaks:
        int_peaks.append((
            mu + np.random.uniform(-0.06, 0.06),
            sig * 1.2,
            amp * np.random.uniform(0.7, 1.2)
        ))
    int_peaks += [
        (11.25,0.022,0.9),
        (7.35,0.020,0.8),
        (6.85,0.020,0.7),
    ]
    intermediate = build_spectrum(ppm, int_peaks)

    # --- APO ---
    apo_peaks = []
    for i, (mu, sig, amp) in enumerate(shared_peaks):
        if i % 5 == 0:
            continue
        apo_peaks.append((
            mu + np.random.uniform(-0.10,0.10),
            sig * 1.5,
            amp * np.random.uniform(0.5,1.0)
        ))
    apo_peaks += [
        (10.15,0.030,0.9),
        (6.75,0.028,0.8),
    ]
    apo = build_spectrum(ppm, apo_peaks)

    return normalize(holo), normalize(intermediate), normalize(apo)


# ==========================================================
# 2D Spectral Model
# ==========================================================

def generate_hsqc_axes(nH=160, nN=120,
                       H_range=(12, 6.5),
                       N_range=(140, 100)):
    """
    Generate HSQC ppm axes.
    """
    H_ppm = np.linspace(*H_range, nH)
    N_ppm = np.linspace(*N_range, nN)
    H, N = np.meshgrid(H_ppm, N_ppm)
    return H_ppm, N_ppm, H, N

def build_2d_spectrum(H, N, peaks):
    S = np.zeros_like(H)
    for muH, muN, sigH, sigN, amp in peaks:
        S += gaussian2D(H, N, muH, muN, sigH, sigN, amp)
    return S

def generate_2d_spectra(H, N, seed=2):

    np.random.seed(seed)

    # Shared HSQC scaffold
    shared_peaks = [
        (11.2,130,0.03,0.6,0.7),
        (10.6,128,0.03,0.6,0.7),
        (9.4,122,0.03,0.6,0.7),
        (9.1,126,0.03,0.6,0.7),

        (8.9,129,0.02,0.8,0.75),
        (8.85,118,0.025,0.6,0.6),
        (8.8,120,0.025,0.5,0.5),
        (8.7,125,0.02,0.5,0.65),
        (8.6,122,0.025,0.48,0.45),
        (8.5,129,0.025,0.7,0.6),
        (8.4,118,0.03,0.6,0.7),
        (8.3,122,0.025,0.6,0.6),
        (8.1,120,0.02,0.8,0.75),

        (8.0,119,0.02,0.55,0.6),
        (7.95,124,0.02,0.75,0.8),
        (7.9,122,0.02,0.8,0.75),
        (7.8,127,0.025,0.6,0.55),
        (7.7,121,0.02,0.65,0.6),
        (7.6,118,0.03,0.5,0.65),
        (7.5,125,0.02,0.75,0.8),

        (7.3,118,0.025,0.7,0.7),
        (7.2,120,0.02,0.75,0.8),
        (7.1,110,0.03,0.6,0.55),
        (6.9,113,0.025,0.65,0.6),
        (6.7,108,0.02,0.7,0.6),
    ]

    # --- HOLO ---
    holo_peaks = shared_peaks + [
        (9.3,128,0.02,0.75,0.8),
        (8.65,135,0.02,0.65,0.55),
        (8.05,127,0.02,0.65,0.55),
        (7.55,112,0.02,0.55,0.6),
    ]
    holo = build_2d_spectrum(H, N, holo_peaks)

    # --- INTERMEDIATE ---
    int_peaks = []
    for muH, muN, sH, sN, amp in shared_peaks:
        int_peaks.append((
            muH + np.random.uniform(-0.06,0.06),
            muN + np.random.uniform(-1.2,1.2),
            sH*1.2,
            sN*1.2,
            amp*np.random.uniform(0.7,1.3)
        ))

    int_peaks += [
        (7.85,129,0.025,0.6,0.6),
        (7.35,119,0.025,0.5,0.6),
        (6.85,106,0.025,0.5,0.7),
    ]

    intermediate = build_2d_spectrum(H, N, int_peaks)

    # --- APO ---
    apo_peaks = []
    for i,(muH,muN,sH,sN,amp) in enumerate(shared_peaks):
        if i % 4 == 0:
            continue
        apo_peaks.append((
            muH + np.random.uniform(-0.1,0.1),
            muN + np.random.uniform(-2,2),
            sH*1.5,
            sN*1.5,
            amp*np.random.uniform(0.5,1.0)
        ))

    apo_peaks += [
        (10.1,121,0.04,0.7,0.65),
        (8.1,125,0.04,0.7,0.65),
        (7.55,126,0.025,0.6,0.6),
        (6.75,117,0.03,0.6,0.6),
    ]

    apo = build_2d_spectrum(H, N, apo_peaks)

    return (normalize(holo),
            normalize(intermediate),
            normalize(apo))
# ==========================================================
# Kinetic Model
# ==========================================================

def populations(t, k1, k2):
    """
    Sequential irreversible:
    Holo -> Intermediate -> Apo
    where k1 ≠ k2
    """
    p_holo = np.exp(-k1 * t)

    p_int = (k1 / (k2 - k1)) * (
        np.exp(-k1 * t) - np.exp(-k2 * t)
    )

    p_apo = 1 - (
        (k2 * np.exp(-k1 * t) - k1 * np.exp(-k2 * t))
        / (k2 - k1)
    )

    return p_holo, p_int, p_apo


def generate_time_vector(n_points=20, dt=15):
    return np.linspace(0, dt * (n_points - 1), n_points)


# ==========================================================
# Mixture Generation
# ==========================================================

def generate_mixture_matrix(holo, intermediate, apo,
                            time_points,
                            k1=0.035,
                            k2=0.011,
                            noise_level=0.01,
                            seed=0):
    """
    Returns:
        D  -> data matrix (time x ppm)
        C_true -> true concentration matrix (time x 3)
    """
    np.random.seed(seed)

    D = []
    C = []

    for t in time_points:
        p_holo, p_int, p_apo = populations(t, k1, k2)

        spectrum = (
            p_holo * holo +
            p_int  * intermediate +
            p_apo  * apo
        )

        sigma = noise_level * np.max(np.abs(spectrum))
        noise = np.random.normal(0, sigma, spectrum.shape)

        D.append(spectrum + noise)
        C.append([p_holo, p_int, p_apo])

    return np.array(D), np.array(C)

def save_dataset(filename, **arrays):
    """
    Save arrays to a compressed NPZ file.

    Example:
    save_dataset(
        "dataset_1D.npz",
        D=D,
        ppm=ppm,
        time_points=time_points
    )
    """
    np.savez_compressed(filename, **arrays)