import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import expit

EXPERIMENTAL_DATA_PATH = os.path.join("input_experimental_files", "input_iv_file.txt")
OUTPUT_DATA_FOLDER = "output_files/"

def load_experimental_data(filepath):
    """Loads experimental data."""
    data_points = []
    try:
        with open(filepath, 'r') as f:
            # Skip the header line
            next(f)
            
            for line in f:
                # 1. strip() removes newline characters
                # 2. split('\t') separates by tab
                # 3. [x for x in ... if x] filters out empty strings (caused by trailing tabs)
                parts = [x for x in line.strip().split('\t') if x]
                
                # Ensure the line has at least 2 columns (Vgs and Ids)
                if len(parts) >= 2:
                    try:
                        # Take the first two columns (Vgs, Ids)
                        v = float(parts[0])
                        i = float(parts[1])
                        data_points.append([v, i])
                    except ValueError:
                        continue # Skip lines that can't be converted to numbers
        
        # Convert to numpy array
        data = np.array(data_points)
        return data[:, 0], data[:, 1]
    
    except FileNotFoundError:
        print(f"Error: Could not find file at {filepath}")
        return np.array([]), np.array([])
    except Exception as e:
        print(f"Error loading data: {e}")
        return np.array([]), np.array([])
    
def calculate_hgte_bandgap(Eg_300K, target_T_K):
    """
    Estimates the bandgap of HgTe NCs at a target temperature using the Varshni relation.
    
    Parameters:
    - Eg_300K (float): The bandgap energy at room temperature (300 K) in eV.
    - target_T_K (float): The target temperature in Kelvin.
    
    Returns:
    - float: The estimated bandgap energy at the target temperature in eV.
    """
    # Varshni parameters for pure HgTe based on Laurenti et al. (1990)
    # alpha is negative because the HgTe bandgap increases with temperature.
    alpha = -6.3e-4  # eV/K 
    beta = 11.0      # K
    
    # 1. Back-calculate the bandgap at 0 K using the known value at 300 K
    # Standard formula: Eg(T) = Eg(0) - (alpha * T^2) / (T + beta)
    # Rearranged: Eg(0) = Eg(T) + (alpha * T^2) / (T + beta)
    Eg_0K = Eg_300K + (alpha * (300**2)) / (300 + beta)
    
    # 2. Calculate the bandgap at the target temperature (e.g., 77 K)
    Eg_target = Eg_0K - (alpha * (target_T_K**2)) / (target_T_K + beta)
    
    return Eg_target


# --- Constants & Simulation Parameters ---
NUM_LAYERS = 2          # Number of layers in the stack
NUM_CURRENT_LAYERS = 2  # Active current-carrying layers at the bottom
SIM_STEPS = 1000       # Sweep resolution

# --- Plotting Window Parameters ---
VGATE_PLOT_MIN = -50  
VGATE_PLOT_MAX = 50 

# --- Disorder Parameters (With Arrhenius Thermal Activation) ---
SIGMA_BROADENING_1Sh = 0.026  # Gaussian energy disorder width (eV)
SIGMA_BROADENING_1Se = 0.047  # Gaussian energy disorder width (eV)
SIGMA_BROADENING_1Pe = 0.072  # Gaussian energy disorder width (eV)
ARRAY_LENGTH = 50000       # Array size for disorder statistics
ALPHA_THERMAL = 1      # Scales the effective activation barrier height (higher = lower barrier)  

# Physics Parameters
N_SH = 10 # hole degeneracy
N_SE = 2  # 1s state degeneracy
N_PE = 6  # 1p state degeneracy

M_SH = 0.01 # 1S_h mobility ((cm^2)/(V*s))
M_SE = 0.048 # 1S_e mobility ((cm^2)/(V*s))
M_PE = 0.056 # 1P_e mobility ((cm^2)/(V*s))

# Energy and Temperature
e = 1.6e-19 # Fundamental charge (C)
E_K = e / 1.38e-23 # e/K_b
T_K = 77 # K
BETA = E_K / T_K           
K_B_T = 1.0 / BETA # Thermal energy (eV)

# Device Dimensions
V_sd = 0.5 
channel_length = 10e-4 
channel_width = 450e-4 
N_electrodes = 9
E_sd = V_sd / channel_length 
NC_radius = 4.8e-7
print(f"Nanocrystal Radius: {NC_radius} cm")
Ligands_length = 0.4e-7 
total_particle_radius = NC_radius + Ligands_length 
total_particle_volume = (4/3)*np.pi*total_particle_radius**3
Q_P = 10.25 # 2.5 charges per 10 dots
E_D = 3.9 # SiO_2 dielectric constant
E_1 = 16.9 # QD film dielectric constant
E_0 = 8.82e-14 #
D_L = 2*total_particle_radius*(np.sqrt(2/3)) 
D_D = 3e-5 
CP_AREA = 2*np.sqrt(3)*total_particle_radius**2

bandgap = calculate_hgte_bandgap(0.31, T_K)
print(f"Calculated bandgap at {T_K}K is {bandgap} eV")
E_1s_center = bandgap/2
E_1p_center = bandgap/2 + 0.20
E_val_center = -bandgap/2

# --- Optimized Multipliers ---
VOLT_MULTIPLIER = (1.6e-19 * D_L)/ (E_1 * E_0 * CP_AREA )
CURR_MULTIPLIER = (e*V_sd*channel_width*(N_electrodes-1))/(2*np.sqrt(3)*(total_particle_radius**2)*channel_length)

# --- Initialize Energy Arrays with Disorder Broadening ---
np.random.seed(42) # Set seed for reproducible disorder profiles
E_val_array = np.random.normal(loc=E_val_center, scale=SIGMA_BROADENING_1Sh, size=ARRAY_LENGTH)
E_1s_array = np.random.normal(loc=E_1s_center, scale=SIGMA_BROADENING_1Se, size=ARRAY_LENGTH)
E_1p_array = np.random.normal(loc=E_1p_center, scale=SIGMA_BROADENING_1Pe, size=ARRAY_LENGTH)

# --- Asymmetric Arrhenius-Modulated Mobility Arrays ---
# Electrons (1S, 1P): E < E_center requires thermal activation (penalty). E >= E_center is spontaneous (max mobility).
M_SE_array = np.where(E_1s_array < E_1s_center, 
                      M_SE * np.exp(-(E_1s_center - E_1s_array) / (ALPHA_THERMAL * K_B_T)), 
                      M_SE)

M_PE_array = M_PE 

# Holes (Valence Band): Trap states for holes are physically located at higher electron energies.
# Therefore, E > E_center requires thermal activation (penalty) for holes. E <= E_center is spontaneous.
M_SH_array = np.where(E_val_array > E_val_center, 
                      M_SH * np.exp(-(E_val_array - E_val_center) / (ALPHA_THERMAL * K_B_T)), 
                      M_SH)


def calculate_layer(fermi_energy_in, voltage_in, E_val_arr, E_1s_arr, E_1p_arr, M_SH_arr, M_SE_arr, M_PE_arr):
    ef_local = fermi_energy_in - voltage_in
    
    s1h = N_SH * expit(-(E_val_arr - ef_local) * BETA)
    s1e = N_SE * expit(-(E_1s_arr - ef_local) * BETA)
    p1e = N_PE * expit(-(E_1p_arr - ef_local) * BETA)
    
    averaged_charge = np.mean(-(s1h + s1e + p1e) + Q_P)
    voltage_out = voltage_in + averaged_charge * VOLT_MULTIPLIER

    # Transport calculation integrated with localized Arrhenius state mobilities
    current_out_arr = ((s1h * (N_SH - s1h) / N_SH) * N_SH * M_SH_arr + 
                       (s1e * (N_SE - s1e) / N_SE) * N_SE * M_SE_arr + 
                       (p1e * (N_PE - p1e) / N_PE) * N_PE * M_PE_arr) * CURR_MULTIPLIER
                       
    averaged_current = np.mean(current_out_arr)
    
    return averaged_current, ef_local, voltage_out, averaged_charge

# =============================================================================
# CALIBRATING RANGE OF FERMI LEVEL VALUES
# =============================================================================

def get_vgate_for_ef(ef_val):
    """Helper function to find the resulting Vgate for a single given initial Ef."""
    current_ef = ef_val
    current_v = 0.0
    for i in range(NUM_LAYERS):
        _, current_ef, current_v, _ = calculate_layer(
            current_ef, current_v, E_val_array, E_1s_array, E_1p_array, M_SH_array, M_SE_array, M_PE_array
        )
    return -current_v * D_D / D_L * E_1 / E_D

def find_ef_boundary(target_vgate):
    """Uses a dynamic binary search to find the exact Ef needed for a target Vgate."""
    low_ef, high_ef = -50.0, 50.0  
    
    for _ in range(10):
        v_low = get_vgate_for_ef(low_ef)
        v_high = get_vgate_for_ef(high_ef)
        if v_low < target_vgate < v_high:
            break
        if v_low >= target_vgate: low_ef -= 20.0
        if v_high <= target_vgate: high_ef += 20.0

    for _ in range(40):
        mid_ef = (low_ef + high_ef) / 2.0
        v_mid = get_vgate_for_ef(mid_ef)
        if v_mid < target_vgate:
            low_ef = mid_ef
        else:
            high_ef = mid_ef
    return (low_ef + high_ef) / 2.0

print("Calibrating simulation boundaries...")
ef_start = find_ef_boundary(VGATE_PLOT_MIN)
ef_end = find_ef_boundary(VGATE_PLOT_MAX)

buffer = (ef_end - ef_start) * 0.05
EF_START_OPTIMIZED = ef_start - buffer
EF_END_OPTIMIZED = ef_end + buffer

# =============================================================================
# SIMULATION RUN
# =============================================================================

print(f"Calibration successful! Sweeping optimized Ef window: [{EF_START_OPTIMIZED:.3f} eV, {EF_END_OPTIMIZED:.3f} eV]")

currents = np.zeros((SIM_STEPS, NUM_LAYERS))
ef_locals = np.zeros((SIM_STEPS, NUM_LAYERS))
voltages = np.zeros((SIM_STEPS, NUM_LAYERS))
charges = np.zeros((SIM_STEPS, NUM_LAYERS))
v_gates = np.zeros(SIM_STEPS)

ef_sweep = np.linspace(EF_START_OPTIMIZED, EF_END_OPTIMIZED, SIM_STEPS)

for t in range(SIM_STEPS):
    current_ef = ef_sweep[t]
    current_v = 0.0
    for i in range(NUM_LAYERS):
        current_i, current_ef, current_v, current_q = calculate_layer(
            current_ef, current_v, E_val_array, E_1s_array, E_1p_array, M_SH_array, M_SE_array, M_PE_array
        )
        currents[t, i], ef_locals[t, i], voltages[t, i], charges[t, i] = current_i, current_ef, current_v, current_q
    
    v_gates[t] = -current_v * D_D / D_L * E_1 / E_D

# --- Masking for Plotting Window ---
mask = (v_gates >= VGATE_PLOT_MIN) & (v_gates <= VGATE_PLOT_MAX)

# Apply mask to all simulation results
v_gates_plot = v_gates[mask]
currents_plot = currents[mask, :]
voltages_plot = voltages[mask, :]
charges_plot = charges[mask, :]
ef_locals_plot = ef_locals[mask, :]

# =============================================================================
# PLOTTING RESULTS
# =============================================================================

# Cosmetics
cmap = plt.get_cmap('viridis')
values = np.linspace(0, 1, NUM_LAYERS)
colors = [cmap(val) for val in values]

# Load experimental data
exp_Vgs, exp_current = load_experimental_data(EXPERIMENTAL_DATA_PATH)

# Plot 1: Charge per dot vs Vgate
plt.figure(figsize=(6, 4))
for i in range(NUM_LAYERS, 0, -1):
    # Access array index directly (Layer 1 is index 0)
    plt.plot(v_gates_plot, charges_plot[:, i-1], linestyle='-', color=colors[NUM_LAYERS - i], label=f'Layer {i}')
plt.xlabel("Vgate (V)")
plt.ylabel("Charge / Dot")
plt.title(f"Layer Charge vs Gate Voltage")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left') 
plt.grid(True)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DATA_FOLDER}charge_per_dot_vs_Vgate.pdf")

# Plot 3: Current vs Vgate (Log Scale)
plt.figure(figsize=(6, 4))
for i in range(NUM_LAYERS, 0, -1):
    plt.plot(v_gates_plot, currents_plot[:, i-1]*1e6, linestyle='-', color=colors[NUM_LAYERS - i], label=f'Layer {i}')
plt.xlabel("Vgate (V)")
plt.ylabel("Current ($\mu$A)")
plt.title(f"Layer Current vs Gate Voltage")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.yscale('log')
plt.ylim(top=1e1, bottom=1e-9)
plt.grid(True)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DATA_FOLDER}layer_current.pdf")

# Plot: Comparison with Experimental Data
# Calculate Total Current
start_layer = NUM_LAYERS - NUM_CURRENT_LAYERS + 1
total_current_plot = np.sum(currents_plot[:, (start_layer-1):NUM_LAYERS], axis=1) * 1e6

# Prepare secondary axis helper
A = N_electrodes * channel_width * channel_length
C_ox = (A * E_D * E_0) / D_D
nc_number = A / (2 * np.sqrt(3) * total_particle_radius**2)

def vgs_to_charge(vgs):
    q_total = (vgs + 7) * C_ox
    charge_number_per_nc = (q_total / nc_number) / e
    return charge_number_per_nc

def charge_to_vgs(charge):
    vgs = (charge * nc_number * e) / C_ox - 7
    return vgs

# Linear Plot
fig1, ax1 = plt.subplots(figsize=(7, 5))
if len(exp_Vgs) > 0:
    ax1.plot(exp_Vgs, exp_current*1e6, '--', label='Experimental', color='green', markersize=3)
ax1.plot(v_gates_plot, total_current_plot, '-', label='Simulated', color='purple', linewidth=1.5)
ax1.set_xlabel('Gate Voltage, $V_{gs}$ (V)')
ax1.set_ylabel('Drain-Source Current, $I_{ds}$ ($\mu$A)')
ax1.grid(True, alpha=0.5)
ax1.legend()
secax1 = ax1.secondary_xaxis('top', functions=(vgs_to_charge, charge_to_vgs))
secax1.set_xlabel('Number of charges per dot')
plt.tight_layout()
plt.savefig(f'{OUTPUT_DATA_FOLDER}plot_linear.pdf')

# Log Plot
fig2, ax2 = plt.subplots(figsize=(7, 5))
if len(exp_Vgs) > 0:
    ax2.plot(exp_Vgs, exp_current*1e6, '--', label='Experimental', color='green', markersize=3)
ax2.plot(v_gates_plot, total_current_plot, '-', label='Simulated', color='purple', linewidth=1.5)
ax2.set_yscale('log')
ax2.set_xlabel('Gate Voltage, $V_{gs}$ (V)')
ax2.set_ylabel('Drain-Source Current, $I_{ds}$ (A)')
ax2.grid(True, which="both", alpha=0.5)
ax2.legend()
secax2 = ax2.secondary_xaxis('top', functions=(vgs_to_charge, charge_to_vgs))
secax2.set_xlabel('Number of charges per dot')
plt.tight_layout()
plt.savefig(f'{OUTPUT_DATA_FOLDER}plot_log.pdf')