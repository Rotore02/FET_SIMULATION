# FET_SIMULATION
This reporitory provides the code for the HgTe QD - based FET simulation.

In what follows we will provide a step-by-step breakdown of the 1D self-consistent electrostatic and hopping transport model implemented to simulate HgTe quantum dot (QD) field-effect transistors (FETs). Section 1 describes the mathematical models and equations used in the simulation. Section 2 describes the code parameters and functions in detail. Section 3 presents a tutorial on how to execute the code.

---

## Section 1: Mathematical Framework & Models

The simulation uses a self-consistent, 1D mean-field electrostatic approximation paired with an ensemble-averaged hopping transport model over a disordered energy landscape. The active layer channel is spatially discretized along the axis perpendicular to the gate dielectric substrate ($z$-axis) into a finite sequence of coupled nanocrystal sheets ($i = 1, 2, \dots, N_{\text{layers}}$).

### 1.1 Structural Geometry & Hexagonal Close-Packed Discrete Sheets
The model assumes that individual colloidal quantum dots are rigid spheres organized in a three-dimensional Hexagonal Close-Packed (HCP) matrix. Let $R_{\text{qd}}$ denote the inorganic nanocrystal core radius, and $R_{\text{lig}}$ represent the insulating ligand shell length. The total interaction radius separating the centers of two adjacent close-packed QDs is given by:

$$R_{\text{tot}} = R_{\text{qd}} + R_{\text{lig}}$$

In an HCP lattice arrangement, the vertical inter-planar layer spacing (thickness of each discrete layer sheet along the $z$-axis) is defined as:

$$D_l = \sqrt{\frac{8}{3}} R_{\text{tot}}$$

The corresponding projected surface footprint area occupied by a single quantum dot within the close-packed horizontal plane is:

$$A_l = 2 \sqrt{3} R_{\text{tot}}^2$$

### 1.2 1D Mean-Field Poisson-Poisson Layer Coupling
Each discrete layer $i$ is treated as a continuous, thin macro-sheet containing a uniform distribution of background doping charges. The average net net charge residing on an individual quantum dot in layer $i$ (expressed in units of elementary charge $e$) is determined by balancing the invariant doping background ($Q_P$) against the mobile electron and hole occupancies across the valence band ($\text{val}$), $1S_e$, and $1P_e$ states:

$$\langle Q^{(i)} \rangle = Q_P - \sum_{k \in \{\text{val}, 1S_e, 1P_e\}} n_k^{(i)}$$

The macroscopic sheet charge density $\sigma_{\text{sheet}}^{(i)}$ characterizing layer $i$ is defined by distributing this net charge across the close-packed cross-sectional footprint area $A_l$:

$$\sigma_{\text{sheet}}^{(i)} = \frac{e \langle Q^{(i)} \rangle}{A_l} = \frac{e \langle Q^{(i)} \rangle}{2 \sqrt{3} R_{\text{tot}}^2}$$

Applying Gauss's Law across a 1D cross-section, the step-wise change in the local electric field $\Delta E^{(i)}$ induced across layer $i$ inside the quantum dot solid (characterized by an effective relative permittivity $\epsilon_1$) is:

$$\Delta E^{(i)} = \frac{\sigma_{\text{sheet}}^{(i)}}{\epsilon_1 \epsilon_0}$$

Integrating this field step across the finite layer thickness $D_l$ yields the incremental change in the cumulative electrostatic potential $V^{(i)}$ from layer $i-1$ to layer $i$:

$$V^{(i)} - V^{(i-1)} = \Delta E^{(i)} D_l = \frac{\sigma_{\text{sheet}}^{(i)}}{\epsilon_1 \epsilon_0} D_l$$

Substituting the geometric expressions for $\sigma_{\text{sheet}}^{(i)}$ and $D_l$ explicitly into this potential step yields:

$$V^{(i)} - V^{(i-1)} = \frac{e \langle Q^{(i)} \rangle}{\epsilon_1 \epsilon_0 \left(2 \sqrt{3} R_{\text{tot}}^2\right)} \cdot \left(\sqrt{\frac{8}{3}} R_{\text{tot}}\right)$$

Isolating and simplifying the algebraic prefactors:

$$\frac{\sqrt{\frac{8}{3}}}{2\sqrt{3}} = \frac{\frac{\sqrt{8}}{\sqrt{3}}}{2\sqrt{3}} = \frac{\sqrt{8}}{2 \cdot 3} = \frac{2\sqrt{2}}{6} = \frac{\sqrt{2}}{3}$$

Thus, the exact recursive potential update expression across the coupled 1D layer grid reads:

$$V^{(i)} = V^{(i-1)} + \langle Q^{(i)} \rangle \frac{\sqrt{2}}{3} \frac{e}{\epsilon_1 \epsilon_0 R_{\text{tot}}}$$

This cumulative potential profile acts as a local electrostatic gate, shifting the local electrochemical environment for carriers. The effective local Fermi level $E_{F,\text{local}}^{(i)}$ governing layer $i$ is shifted relative to the initial injection reference Fermi level $E_{F,0}$ at the channel interface:

$$E_{F,\text{local}}^{(i)} = E_{F,0} - e V^{(i-1)}$$

### 1.3 Boundary Conditions & Macroscopic Gate Voltage Mapping
To close the self-consistent loop and relate the boundary potential of the multi-layer film to the physical potential applied at the external gate electrode ($V_G$), the model enforces the continuity of the electric displacement field vector $\mathbf{D}$ across the interface between the last QD layer ($i = N_{\text{layers}}$) and the underlying gate oxide dielectric (defined by thickness $t_{\text{ox}}$ and relative permittivity $\epsilon_d$):

$$D_{\text{ox}} = D_{\text{film}} \implies \epsilon_d \epsilon_0 E_{\text{ox}} = \epsilon_1 \epsilon_0 E_{\text{film}}$$

The electric field across the gate oxide is driven by the potential difference between the terminal gate contact and the last layer of the film. By rewriting the previous relation as a function of the local voltage $V^{(N_{layers})}$ and the gate voltage $V_G$, one obtains:

$$V_G = \frac{\epsilon_1}{\epsilon_d} \sqrt{\frac{3}{8}} \frac{t_{\text{ox}}}{R_{\text{tot}}} V^{(N_{\text{layers}})}$$

which is the $V_G$ that needs to be applied to obtain the considered charge and potential distribution across the film.

### 1.4 Energetic Disorder Landscapes & Density of States (DOS)
To capture structural non-idealities, such as quantum dot size polydispersity, the simulation maps an energetic disorder profile. The density of states for each quantum level is represented by a discrete ensemble of $M$ independent, non-interacting sub-vessels or sites. For each state $k \in \{\text{val}, 1S_e, 1P_e\}$, the localized site energy levels $E_{k,j}$ (where $j = 1, 2, \dots, M$) are sampled from a normal Gaussian distribution centered around the nominal energy level $E_k^0$ with a broadening standard deviation $\sigma_k$:

$$P(E_{k,j}) = \frac{1}{\sqrt{2\pi}\sigma_k} \exp\left( -\frac{(E_{k,j} - E_k^0)^2}{2\sigma_k^2} \right)$$

### 1.5 Fermi-Dirac Distribution and Thermal Carrier Statistics
Because individual layers are assumed to maintain local thermodynamic equilibrium, the probability $f$ that a carrier populates a specific disordered energy site $j$ belonging to state $k$ within layer $i$ is strictly governed by Fermi-Dirac statistics:

$$f(E_{k,j}, E_{F,\text{local}}^{(i)}) = \frac{1}{1 + \exp\left( \frac{E_{k,j} - E_{F,\text{local}}^{(i)}}{k_B T} \right)}$$

The total number of accumulated mobile carriers $n_k^{(i)}$ in state $k$ inside layer $i$ is calculated by integrating across the disordered landscape, which is evaluated numerically by taking the ensemble mean over the $M$ sites and multiplying by the structural state degeneracy factor $N_k$:

$$n_k^{(i)} = \frac{N_k}{M} \sum_{j=1}^{M} f(E_{k,j}, E_{F,\text{local}}^{(i)})$$

### 1.6 Macroscopic Charge Transport
In a quantum-confined nanocrystal array, transport occurs via a series of sequential, phonon-assisted quantum mechanical tunneling steps (hopping) between adjacent localized dot sites. Under the Pauli Exclusion Principle, a charge carrier can successfully hop to a target dot site only if that state is vacant. Consequently, the net hopping probability is proportional to the product of the occupied site probability and the vacancy availability factor, yielding an expression scaled by $n_{k,j}^{(i)} \left(N_k - n_{k,j}^{(i)}\right)$, where $n_{k,j}^{(i)}$ is the amount of charges in quantum state $k$, disorder site $j$ and layer $i$, while $N_k$ is the multiplicity of state $k$.

The current density $J_{k,j}^{(i)}$ contributed by the quantum state $k$, disorder site $j$ in layer $i$ is modeled as:

$$J_{k,j}^{(i)} = e (N_{\text{finger}} - 1) \rho_{k,j}^{(i)} v_{k,j}$$

where $N_{\text{finger}}$ is the number of interdigitated finger pairs, and $\rho_{k,j}^{(i)}$ is the localized carrier density expressed as:

$$\rho_{k,j}^{(i)} = \frac{n_{k,j}^{(i)} \left(N_k - n_{k,j}^{(i)}\right)}{A_l}$$

The carrier drift velocity $v_{k,j}$ under a small source-drain bias field $E_{\text{sd}}$ is related through the localized mobility $\mu_{k,j}$ by $v_{k,j} = \mu_{k,j} E_{\text{sd}}$, where $E_{\text{sd}} = \frac{V_{\text{sd}}}{L}$, and $L$ represents the longitudinal channel length. 

Integrating this current density across the vertical cross-sectional sheet area of a single close-packed layer ($A_{\text{sheet}} = W \cdot D_l$, where $W$ is the channel width) yields:

$$I_{k,j}^{(i)} = J_{k,j}^{(i)} A_{\text{sheet}} = e (N_{\text{finger}} - 1) \left[ \frac{n_{k,j}^{(i)} \left(N_k - n_{k,j}^{(i)}\right)}{2\sqrt{3}R_{\text{tot}}^2} \right] \mu_{k,j} \left(\frac{V_{\text{sd}}}{L}\right) \left(W \cdot \sqrt{\frac{8}{3}}R_{\text{tot}}\right)$$

Grouping the structural variables into a single unified global prefactor constant $\alpha_0$:

$$\alpha_0 = \left[ e (N_{\text{finger}} - 1) \frac{W V_{\text{sd}}}{L \cdot R_{\text{tot}}^2} \right] \cdot \frac{\sqrt{\frac{8}{3}}}{2\sqrt{3}} = \frac{\sqrt{3}}{6} \frac{W e (N_{\text{finger}} - 1) V_{\text{sd}}}{R_{\text{tot}}^2 L}$$

The total macroscopic current traversing layer $i$, denoted as $I^{(i)}$, is found by summing over all valid energy levels $k$ and calculating the ensemble average over the sample array of size $M$:

$$I^{(i)} = \frac{\alpha_0}{M} \sum_{j=1}^{M} \sum_{k \in \{\text{val}, 1S_e, 1P_e\}} n_{k,j}^{(i)} \left(N_k - n_{k,j}^{(i)}\right) \mu_{k,j}$$

The net macroscopic drain current ($I_D$) flowing through the active region of the transistor is determined by summing these layer currents over the designated subset of conducting channels ($\mathcal{S}_{\text{current}}$) near the transport plane:

$$I_D = \sum_{i \in \mathcal{S}_{\text{current}}} I^{(i)}$$

### 1.7 Microscopic Miller-Abrahams (Arrhenius) Mobility Modulation
The baseline mobility tracking parameters ($\mu_0^k$) are modulated using a microscopic Miller-Abrahams formulation to account for the thermal activation energy required to escape deep energetic traps:

* **Valence State Band ($k = \text{val}$):** Traps for holes correspond to higher electron energies. Hops to states residing above the nominal band center energy ($E_{\text{val},j} > E_{\text{val}}^0$) face an exponential Arrhenius penalty:
    $$\mu_{\text{val},j} = \begin{cases} \mu_0^{\text{val}} \exp\left( -\frac{E_{\text{val},j} - E_{\text{val}}^0}{\alpha_T k_B T} \right) & \text{for } E_{\text{val},j} > E_{\text{val}}^0 \\ \mu_0^{\text{val}} & \text{for } E_{\text{val},j} \le E_{\text{val}}^0 \end{cases}$$

* **$1S_e$ Conduction State Band ($k = 1S_e$):** Traps for conduction electrons correspond to localized sites located below the nominal band center. Hops from states where $E_{1S_e,j} < E_{1S_e}^0$ are thermally restricted:
    $$\mu_{1S_e,j} = \begin{cases} \mu_0^{1S_e} \exp\left( -\frac{E_{1S_e}^0 - E_{1S_e,j}}{\alpha_T k_B T} \right) & \text{for } E_{1S_e,j} < E_{1S_e}^0 \\ \mu_0^{1S_e} & \text{for } E_{1S_e,j} \geq E_{1S_e}^0 \end{cases}$$

* **$1P_e$ Conduction State Band ($k = 1P_e$):** Modeled as a highly delocalized transport state with minimal structural trapping, maintaining uniform baseline mobility across all configurations:
    $$\mu_{1P_e,j} = \mu_0^{1P_e}$$

---

## Section 2: Code Documentation and Implementation Details

This section details the design, parameter requirements, and structural functions contained within the Python simulation script.

### 2.1 Folder Tree Structure 

```text
FET_SIMULATION
├── input_experimental_files
│   └── input_iv_file.txt
├── output_files
│   ├── charge_per_dot_vs_Vgate.pdf
│   ├── layer_current.pdf
│   ├── plot_linear.pdf
│   └── plot_log.pdf
├── FET_simulation.py
└── README.md
```

### 2.2 Unified Model Parameter Definition
The variables governing the execution of the simulation script are detailed below, categorized by their physical domains.

| Variable Name | Mathematical Symbol | Description | Physical Units |
| :--- | :--- | :--- | :--- |
| **Simulation Setup** | | | |
| `NUM_LAYERS` | $N_{\text{layers}}$ | Total number of discrete quantum dot layers within the film | Dimensionless |
| `NUM_CURRENT_LAYERS` | $\mathcal{S}_{\text{current}}$ | Subset of active layers participating in lateral current transport | Dimensionless |
| `SIM_STEPS` | — | Total resolution points computed during the Fermi energy loop sweep | Dimensionless |
| `VGATE_PLOT_MIN` | $V_{G,\text{min}}$ | Lower boundary limit of the gate voltage plotting window | V |
| `VGATE_PLOT_MAX` | $V_{G,\text{max}}$ | Upper boundary limit of the gate voltage plotting window | V |
| `ARRAY_LENGTH` | $M$ | Sampling scale size used to construct the Gaussian disorder arrays | Dimensionless |
| **Disorder Parameters** | | | |
| `SIGMA_BROADENING_1Sh` | $\sigma_{\text{val}}$ | Standard deviation parameter for Gaussian valence band disorder | eV |
| `SIGMA_BROADENING_1Se` | $\sigma_{1S_e}$ | Standard deviation parameter for Gaussian $1S_e$ conduction disorder | eV |
| `SIGMA_BROADENING_1Pe` | $\sigma_{1P_e}$ | Standard deviation parameter for Gaussian $1P_e$ conduction disorder | eV |
| `ALPHA_THERMAL` | $\alpha_T$ | Empirical scaling factor for Miller-Abrahams hopping activation | Dimensionless |
| **Degeneracy & Mobility** | | | |
| `N_SH` | $N_{\text{val}}$ | Electronic state degeneracy factor for the valence level | Dimensionless |
| `N_SE` | $N_{1S_e}$ | Electronic state degeneracy factor for the $1S_e$ conduction level | Dimensionless |
| `N_PE` | $N_{1P_e}$ | Electronic state degeneracy factor for the $1P_e$ conduction level | Dimensionless |
| `M_SH` | $\mu_0^{\text{val}}$ | Intrinsic baseline low-field mobility of the valence state | $\text{cm}^2 \cdot \text{V}^{-1} \cdot \text{s}^{-1}$ |
| `M_SE` | $\mu_0^{1S_e}$ | Intrinsic baseline low-field mobility of the $1S_e$ state | $\text{cm}^2 \cdot \text{V}^{-1} \cdot \text{s}^{-1}$ |
| `M_PE` | $\mu_0^{1P_e}$ | Intrinsic baseline low-field mobility of the $1P_e$ state | $\text{cm}^2 \cdot \text{V}^{-1} \cdot \text{s}^{-1}$ |
| **Device & Physics** | | | |
| `e` | $e$ | Fundamental elementary electron charge constant | C |
| `T_K` | $T$ | Simulation environment absolute operational temperature | K |
| `V_sd` | $V_{\text{DS}}$ | Fixed applied drain-to-source bias potential | V |
| `channel_length` | $L$ | Total physical path length of the active channel | cm |
| `channel_width` | $W$ | Transverse spatial width of the active channel | cm |
| `N_electrodes` | $N_{\text{finger}}$ | Number of fingers in the interdigitated electrode structure | Dimensionless |
| `NC_radius` | $R_{\text{qd}}$ | Inorganic radius of the nanocrystal core | cm |
| `Ligands_length` | $R_{\text{lig}}$ | Longitudinal length of the capping organic ligand chain | cm |
| `Q_P` | $Q_P$ | Background core structural ionized chemical doping charge | Charges / dot |
| `E_D` | $\epsilon_d$ | Relative permittivity (dielectric constant) of the oxide substrate | Dimensionless |
| `E_1` | $\epsilon_1$ | Relative effective permittivity of the assembled QD solid film | Dimensionless |
| `t_ox` | $t_{\text{ox}}$ | Thickness profile of the underlying gate oxide substrate layer | cm |

### 2.3 Functional Architecture & Core Methods Breakdown

#### 1. `load_experimental_data(filepath)`
* **Logical Operations:**
  1. Opens the target data file via an institutional file-stream reader, skipping the initial header text index line.
  2. Sequentially parses rows, stripping line breaks, and splits characters using standard tab delimiters (`\t`).
  3. Trailing white space and null elements are scrubbed via a list comprehension filter.
  4. Segregates the array data matrix: Column 0 maps to the experimental Gate Voltage vector ($V_{gs}$), and Column 1 captures the experimental Drain Current data ($I_{ds}$).
* **Returns:** A nested tuple containing two independent NumPy float arrays: `(vgate_experimental, idrain_experimental)`.

#### 2. `calculate_hgte_bandgap(Eg_300K, target_T_K)`
* **Logical Operations:**
  1. Standardizes empirical Varshni material constants calibrated for bulk mercury telluride systems ($\alpha_{\text{Varshni}} = 0.00033\text{ eV/K}$, $\beta_{\text{Varshni}} = 160.0\text{ K}$).
  2. Back-calculates the theoretical material bandgap boundary profile normalized to absolute zero ($E_g(0\text{ K})$) from the user-provided room-temperature index entry (`Eg_300K`):
     $$E_g(0) = E_g(300\text{ K}) - \frac{\alpha_{\text{Varshni}} \cdot 300^2}{300 + \beta_{\text{Varshni}}}$$
  3. Projects and computes the expected baseline bandgap configuration matching the operational temperature runtime parameter (`target_T_K`):
     $$E_g(T) = E_g(0) + \frac{\alpha_{\text{Varshni}} T^2}{T + \beta_{\text{Varshni}}}$$
* **Returns:** A floating-point number representing the temperature-adjusted energy bandgap value (eV).

#### 3. `calculate_layer(fermi_energy_in, voltage_in, E_val_arr, E_1s_arr, E_1p_arr, M_SH_arr, M_SE_arr, M_PE_arr)`
* **Logical Operations:**
  1. Computes the layer-specific local Fermi level by tracking the potential drop across the spatial boundary: `ef_local = fermi_energy_in - voltage_in`.
  2. Evaluates the site-by-site thermodynamic state occupancy vector over the entire length $M$ of the disorder array. To optimize execution speed, it uses the high-performance logistic sigmoid mapping function `scipy.special.expit()`, which evaluates the Fermi-Dirac integral:
     $$f = \text{expit}\left(-\frac{E_{j} - E_{F,\text{local}}}{k_B T}\right) = \frac{1}{1 + \exp\left(\frac{E_{j} - E_{F,\text{local}}}{k_B T}\right)}$$
  3. Computes the ensemble average net charge per dot $\langle Q^{(i)} \rangle$ across the $M$ sites, combining the structural ionized background core charge with the mobile electronic state occupancies: `averaged_charge = Q_P - (mean(s1h_occupancies) + mean(s1e_occupancies) + mean(p1e_occupancies))`.
  4. Updates the cumulative electrostatic potential to be passed to the next layer in the loop sequence using the 1D discrete Poisson equation derivation: `voltage_out = voltage_in + averaged_charge * VOLT_MULTIPLIER`.
  5. Evaluates the hopping current across the $M$ configuration paths using the state congestion relation $n(N-n)$, scaled by the local mobility arrays and the unified geometric prefactor:
     $$\text{current\_array} = \alpha_0 \cdot \sum n_k \cdot (N_k - n_k) \cdot \mu_k$$
  6. Averages the array values across the disorder landscape to extract the single macroscopic layer current contribution: `averaged_current = np.mean(current_array)`.
* **Returns:** A tuple of computed layer values: `(averaged_current, ef_local, voltage_out, averaged_charge)`.

#### 4. `get_vgate_for_ef(ef_val)`
* **Logical Operations:**
  1. Initializes a sequential layer-by-layer loop execution block starting from zero initial potential boundary metrics.
  2. Iterates through all $N_{\text{layers}}$ layers using the `calculate_layer` routine for a fixed guess of the Fermi Level ($E_{F,0}$).
  3. Captures the potential profile generated at the final oxide interface layer, $V^{(N_{\text{layers}})}$.
  4. Projects this boundary condition to output the corresponding external terminal macroscopic Gate Voltage ($V_G$) using the continuity of the displacement field.
* **Returns:** A floating-point value representing the macroscopic gate voltage $V_G$ matching the input $E_{F,0}$ level.

#### 5. `find_ef_boundary(target_vgate)`
* **Logical Operations:**
  1. Implements a binary search routine to find the Fermi Level ($E_{F,0}$) required to produce a specific target gate voltage.
  2. *Phase 1 (Dynamic Bracketing):* Initializes search limits and expands the bounds (`low_ef`, `high_ef`) outward until the target $V_G$ is successfully bracketed.
  3. *Phase 2 (Bisection Optimization):* Executes 40 refinement iterations. In each step, it calculates the midpoint `mid_ef = (low_ef + high_ef) / 2`, calls `get_vgate_for_ef(mid_ef)`, and adjusts the search window bounds accordingly.
* **Returns:** A floating-point value specifying the optimized Fermi Level $E_{F,0}$ that yields `target_vgate`.

### 2.4 Execution Sequence & Main Simulation Loop Flow
1. **Disorder Grid Generation:** The script samples from normal Gaussian distributions to generate static energy arrays of size $M$ for the `E_val`, `E_1S_e`, and `E_1P_e` bands, mapping the energetic disorder landscape.
2. **Miller-Abrahams Mobility Initialization:** Constructs static mobility vectors for each site. It applies an asymmetric Arrhenius penalty to the mobilities based on each site's energy position relative to its nominal band center ($E_k^0$).
3. **The Discrete Fermi Sweep Loop:** The program sets up a linear loop split into `SIM_STEPS` intervals across the calibrated energy window. At each index step:
   * It sweeps sequentially through the $N_{\text{layers}}$ coupled film sheets.
   * The cumulative output potential `voltage_out` calculated from layer $i-1$ is piped as the baseline input potential `voltage_in` for layer $i$.
4. **Macroscopic Property Extraction:** Aggregates the simulation metrics. The total drain current ($I_D$) is calculated by summing the layer currents across the designated active subset layer index array (`NUM_CURRENT_LAYERS`).
5. **Data Visualization & Export:** Applies data masks to filter out any numerical boundary errors. It outputs the final simulation results alongside the loaded experimental curves, producing linear and logarithmic transfer characteristics ($I_D$ vs. $V_G$), layer-resolved current profiles, and charge-injection scaling tracks. The plots are saved in .pdf format in the `output_files` folder.

## Section 3: Tutorial

1. Download (clone) the repository to your local machine by running the following command in your terminal or command prompt:
   ```bash
   git clone [https://github.com/Rotore02/FET_SIMULATION.git](https://github.com/Rotore02/FET_SIMULATION.git)
   ```
   *(Note: After cloning, navigate into the downloaded folder using `cd FET_SIMULATION`)*

2. Put the input experimental data file inside the `input_experimental_files` folder. This file needs to be organized in two columns, the first will be the gate voltage, the second will be the source drain current (note that the first row will be skipped by the file reading function). The name of this input file needs to be updated in the `FET_simulation.py` code, precisely in the `EXPERIMENTAL_DATA_PATH` variable.

3. Select the simulation parameters in the `FET_simulation.py` file.

4. Run the simulation using the following commands in your terminal or command prompt:

   **From inside the `FET_SIMULATION` folder:**

   * **Windows:**
     ```cmd
     python FET_simulation.py
     ```

   * **Linux / macOS:**
     ```bash
     python3 FET_simulation.py
     ```

   **From outside the `FET_SIMULATION` folder:**

   * **Windows:**
     ```cmd
     python FET_SIMULATION\FET_simulation.py
     ```

   * **Linux / macOS:**
     ```bash
     python3 FET_SIMULATION/FET_simulation.py
     ```

5. The plots will be saved in .pdf format in the `output_files` folder.