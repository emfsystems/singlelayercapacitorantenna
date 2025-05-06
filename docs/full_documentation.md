# 📚 Full Documentation: SingleLayerCapacitorAntenna\_new Class

---

## 🔍 Overview

The `SingleLayerCapacitorAntenna_new` class models a **single-layer capacitor-based antenna** structure for RF and wireless system simulations. It provides tools to:

* Construct single elements, rectangular arrays, and compound subarrays
* Model microstrip transmission feed lines
* Calculate input impedance, S-parameters, and VSWR
* Simulate far-field radiation patterns
* Handle ground-plane reflections (using Image Theory)
* Simulate receive-mode operation (incoming wave to feed voltages)

This class is engineered for both **academic research** and **industry applications**.

---

## 🔄 Class Architecture

### Core Properties

| Property                                            | Description                                     |
| :-------------------------------------------------- | :---------------------------------------------- |
| `NumberOfElements`                                  | Total number of antenna elements                |
| `ElementLength`, `ElementWidth`, `ElementThickness` | Physical dimensions of each antenna element     |
| `OperatingFrequency`                                | Default operating frequency (Hz)                |
| `FeedVoltage`, `SourceImpedance`                    | Feeding configuration parameters                |
| `Capacitance`, `Resistance`, `Inductance`           | RLC model for single-layer capacitor            |
| `ArrayLength`, `ArrayWidth`                         | Overall array dimensions                        |
| `GroundPlane`                                       | Boolean toggle for ground reflection modeling   |
| `CompoundArrayDefinitions`                          | Defines subarray structures for compound arrays |
| `MultiPort`                                         | Enables multi-feed excitation for arrays        |

### Core Methods

| Method                                                   | Purpose                                                  |
| :------------------------------------------------------- | :------------------------------------------------------- |
| `show()`                                                 | Visualize the antenna geometry                           |
| `lineBasedImpedance(freq, idx)`                          | Compute input impedance using microstrip model           |
| `calculateVSWR(freq)`                                    | Calculate voltage standing wave ratio                    |
| `calculateReturnLoss(freq)`                              | Calculate return loss                                    |
| `visualizeRadiationPattern()`                            | Plot 2D radiation pattern                                |
| `visualize3DPattern()`                                   | Plot 3D radiation pattern                                |
| `currentDistribution(freq)`                              | Solve for current distribution across elements           |
| `simulateReception(incidentField, freq)`                 | Simulate receive voltage response to incoming plane wave |
| `exportSParameters(filename)`                            | Save S11 or S-parameter data as .s1p files               |
| `estimateBandwidth(varargin)`                            | Estimate usable bandwidth based on VSWR threshold        |
| `performImagingScan(azGrid, elGrid, freq, polarization)` | Perform scanning reception simulation over grid angles   |

---

## 👩‍💻 Usage Workflow

### Basic Setup (Single Element)

```matlab
antenna = SingleLayerCapacitorAntenna_new(1, 0.02, 0.001, 0.03, 0.03, ...
    2.4e9, 1, 470e-12, 4*pi*1e-7, 100, 1e-9, 1, 50, 'ConfigurationType', 'single');

antenna.show();
VSWR = antenna.calculateVSWR(2.4e9);
```

### Rectangular Array

```matlab
antenna = SingleLayerCapacitorAntenna_new(8, 0.02, 0.001, 0.06, 0.04, ...
    2.4e9, 1, 470e-12, 4*pi*1e-7, 100, 1e-9, 1, 50, ...
    'ConfigurationType', 'rectangular', 'NumRows', 2, 'NumCols', 4);
```

### Compound Subarrays

```matlab
% Define four 2x2 subarrays
subArrays(1).NumRows = 2; subArrays(1).NumCols = 2;
subArrays(1).ArrayOffset = [-0.03, -0.03, 0];
... % define remaining subarrays
antenna = SingleLayerCapacitorAntenna_new(16, 0.02, 0.001, 0.04, 0.04, ...
    2.4e9, 1, 470e-12, 4*pi*1e-7, 100, 1e-9, 1, 50, 'ConfigurationType', 'compound', ...
    'CompoundArrayDefinitions', subArrays, 'MultiPort', true);
```

---

## 🔢 Key Equations

### Capacitance Reactance

$$
X_c = \frac{1}{2 \pi f C}
$$

### Effective Input Impedance (Including Feed Line)

$$
Z_{in} = Z_0 \frac{Z_{load} + j Z_0 \tan(\beta l)}{Z_0 + j Z_{load} \tan(\beta l)}
$$

Where:

* \$Z\_0\$ = microstrip characteristic impedance
* \$Z\_{load}\$ = load impedance (capacitor RLC equivalent)
* \$\beta\$ = phase constant
* \$l\$ = feed line length

### Radiation Efficiency

$$
\eta = \frac{R_r}{R_r + R_{loss}}
$$

Where:

* \$R\_r\$ = radiation resistance
* \$R\_{loss}\$ = ESR and microstrip loss contributions

### Receive Mode Voltage (Loaded)

$$
V_{loaded} = \frac{E_{incident} \cdot Effective\_Aperture}{Z_{load}}
$$

---

## 📊 Example Output Figures

* **Radiation Pattern:**

  ![Example Radiation Pattern](../figures/example_compound_pattern.png)

* **S11 Reflection Coefficient (VSWR):**

  *Coming Soon*

---

## 🌐 References

* EMF Systems, "Revolutionizing Wireless Communication with Single Layer Capacitor-Based Antenna Technology," arXiv:2403.01170.
* Pozar, David M., "Microwave Engineering," Wiley, 4th Edition.
* Balanis, Constantine A., "Antenna Theory: Analysis and Design," Wiley, 4th Edition.

For additional insights, visit:

* [EMF Systems Official Website](https://emfsystems.africa)

---

# 🔹 End of Documentation
