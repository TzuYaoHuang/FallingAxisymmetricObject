
"""Utility script for estimating terminal velocities and related quantities for
simple axisymmetric shapes falling in water.

All inputs use CGS units: length in cm, mass in g, time in s.
"""

import numpy as np
from typing import Dict

# === Physical constants (CGS) ===
G: float = 980.665  # cm/s^2
DENSITY_WATER: float = 1.0  # g/cm^3
VISCOSITY_DEFAULT: float = 0.01  # g/(cm·s)


def terminal_velocity(m_excess: float, Cd: float, A: float, density: float = DENSITY_WATER, g: float = G) -> float:
    """Return the terminal velocity (cm/s) for a given excess mass, drag coeff and area.
    """
    return float(np.sqrt(2 * g * m_excess / (density * A * Cd)))


def time_to_fall(distance: float, v_t: float) -> float:
    """Estimate time to fall a given distance (cm) at constant velocity v_t (cm/s)."""
    return distance / v_t


def time_to_terminal(m_excess: float, Cd: float, A: float, density: float = DENSITY_WATER, g: float = G, fraction: float = 0.99) -> float:
    """Time to reach a given fraction of terminal velocity (s).

    Uses the analytical solution for motion with quadratic drag.
    """
    v_t = terminal_velocity(m_excess, Cd, A, density, g)
    return float(v_t / g * np.arctanh(fraction))


def distance_to_terminal(m_excess: float, Cd: float, A: float, density: float = DENSITY_WATER, g: float = G, fraction: float = 0.99) -> float:
    """Distance (cm) required to reach a fraction of terminal velocity."""
    v_t = terminal_velocity(m_excess, Cd, A, density, g)
    t_t = time_to_terminal(m_excess, Cd, A, density, g, fraction)
    return float(v_t**2 / g * np.log(np.cosh(g * t_t / v_t)))


def reynolds_number(velocity: float, diameter: float, density: float = DENSITY_WATER, viscosity: float = VISCOSITY_DEFAULT) -> float:
    """Compute Reynolds number (dimensionless) using CGS units."""
    return float(density * velocity * diameter / viscosity)


# === Default parameters ===
H = 55.0  # fall distance in cm
DIAMETER_DEFAULT = 5.0  # cm (used for volume calculations and area)
A = float(np.pi * (DIAMETER_DEFAULT / 2) ** 2)  # cross-sectional area, cm^2
m_excess = 1.8  # g (excess mass)

# Drag coefficients for different shapes
CD_VALUES: Dict[str, float] = {
    "Sphere": 0.45,
    "Ellipsoid": 0.30,
    "Bullet": 0.65,
    "NACA": 0.15,
}

# Characteristic lengths for Reynolds number calculation (cm)
LENGTHS: Dict[str, float] = {
    "Sphere": 5.0,
    "Ellipsoid": 7.5,
    "Bullet": 7.5,
    "NACA": 16.1,
}

# === Integral constant for NACA thickness distribution ===
I_NACA = 0.005507744279444381  # ∫ f(ξ)^2 dξ for 4-digit NACA thickness


def volume(diameter: float, shape: str, aspect_ratio: float = 1.5, naca_t: float = 0.31) -> float:
    """Return volume (cm^3) for a given axisymmetric shape and max diameter.

    Supported shapes: 'Sphere', 'Ellipsoid', 'Bullet', 'NACA'.
    """
    shape_key = shape.title()
    if shape_key == "Sphere":
        return float(np.pi * diameter ** 3 / 6.0)
    if shape_key == "Ellipsoid":
        return float((np.pi * aspect_ratio / 6.0) * diameter ** 3)
    if shape_key == "Bullet":
        return float(np.pi * diameter ** 3 / 3.0)
    if shape_key == "Naca":
        return float((25.0 * np.pi * I_NACA / naca_t) * diameter ** 3)
    raise ValueError(f"Unknown shape: {shape}")


def main() -> None:
    """Print a short report for default shapes and parameters."""
    print(
        f"Falling {H:.0f} cm with cross-sectional area {A:.2f} cm^2 and excess mass {m_excess} g in water (density {DENSITY_WATER} g/cm^3):\n"
    )

    for shape, cd in CD_VALUES.items():
        print("=" * 80)
        v_t = terminal_velocity(m_excess, cd, A)
        t_fall = time_to_fall(H, v_t)
        print(f"{shape}: time to fall {H:.0f} cm at terminal velocity: {t_fall:.2f} s")

        V = volume(DIAMETER_DEFAULT, shape)
        min_mass = V * DENSITY_WATER
        print(f"Minimum mass to sink for {shape}: {min_mass:.2f} g")

        t_95 = time_to_terminal(m_excess, cd, A, fraction=0.95)
        d_95 = distance_to_terminal(m_excess, cd, A, fraction=0.95)
        t_99 = time_to_terminal(m_excess, cd, A)
        d_99 = distance_to_terminal(m_excess, cd, A)
        print(
            f"( {t_95:.2f} s / {d_95:.2f} cm to reach 95% of v_t  |  {t_99:.2f} s / {d_99:.2f} cm to reach 99% of v_t )"
        )

        Re = reynolds_number(v_t, DIAMETER_DEFAULT)
        print(f"Terminal velocity: {v_t:.2f} cm/s")
        print(f"Reynolds number at terminal velocity: {Re:.0f}\n")


if __name__ == "__main__":
    main()