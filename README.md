# 🛁 SingleLayerCapacitorAntenna\_new Class Repository

**Compact, High-Efficiency Antenna Modeling and Simulation Toolkit**

&#x20;&#x20;

---

## 🚀 Introduction

The `SingleLayerCapacitorAntenna_new` MATLAB class provides an end-to-end framework for modeling, analyzing, and simulating **single-layer capacitor-based antennas (SLC antennas)** — a disruptive technology in wireless communication miniaturization.

With capabilities ranging from microstrip feed modeling to full radiation pattern computation and reception simulation, this library supports a wide range of academic, industrial, and research applications.

This repository is actively maintained by [**EMF Systems**](https://emfsystems.africa).

---

## 🎯 Key Features

* 🛁 **Single Element and Multi-Element Array Modeling**
* 🔗 **Microstrip Feed Line and Lumped Element Integration**
* 🌎 **Ground-Plane Reflection (Image Theory) Support**
* 📈 **S-Parameter and Impedance Analysis**
* 🧐 **Receive-Mode Simulation for Incoming Plane Waves**
* 📊 **3D Radiation Patterns, Directivity, Gain, and Efficiency Computation**
* 📚 **Exportable Pattern Data and S-Parameter Files (.s1p)**
* 🛠️ **Extendable for WiFi, 5G, Satellite, IoT Communication Systems**

---

## 🛠 Installation

* MATLAB **R2021b** or later is recommended
* Optional but recommended: **Antenna Toolbox**
* Clone the repository:

  ```bash
  git clone https://github.com/emfsystems/SingleLayerCapacitorAntenna.git
  ```
* Add the `src/` folder to your MATLAB path.

---

## 🧪 Usage Examples

| Example                                                                                      | Description                                            |
| -------------------------------------------------------------------------------------------- | ------------------------------------------------------ |
| [`examples/basic/example1_single_element.m`](examples/basic/example1_single_element.m)       | Single antenna setup (2.4 GHz) with radiation analysis |
| [`examples/basic/example2_rectangular_array.m`](examples/basic/example2_rectangular_array.m) | 8-element 2×4 array configuration                      |
| [`examples/basic/example3_compound_array.m`](examples/basic/example3_compound_array.m)       | 4 subarrays (2×2 elements each), multi-port feeding    |

> See more upcoming system-level integrations under [`research/`](research/).

---

## 📚 Full Documentation

Detailed class documentation is available here:

* [`docs/Full_Documentation.md`](docs/Full_Documentation.md)

Covers:

* Class architecture
* Object properties and methods
* Usage workflows
* Equations and physical modeling
* Code snippets

---

## 🌍 Technology Impact

The SLC Antenna has been shown to deliver:

* **High radiation efficiency (97–99%)**
* **Wide bandwidth (kHz–20 GHz)**
* **Miniaturization advantages for wireless modules**

See the market and technological impact analysis:

* [`docs/Impact_Analysis_Summary.md`](docs/Impact_Analysis_Summary.md)

---

## 🔥 Upcoming Extensions

This repository is growing into a major resource for SLC antenna integration into real-world wireless technologies:

* WiFi access points and mesh networks
* 5G massive MIMO systems
* Satellite communications (LEO, MEO, GEO)
* IoT gateway devices
* LPWAN (LoRa, NB-IoT) systems
* radio astronomy and deep space communication systems

View the development roadmap:

* [`research/planned_studies.md`](research/planned_studies.md)

---

## 🏷 License

Distributed under the **MIT License**.
See [`LICENSE`](LICENSE) for more information.

---

## 🧠 Citation

If you use this library in your research, please cite:

```
@misc{slc_antenna_2025,
  author = {EMF Systems},
  title = {Single Layer Capacitor-Based Antenna Modeling Toolkit},
  year = {2025},
  url = {https://github.com/emfsystems/SingleLayerCapacitorAntenna}
}
```

You can also reference our white paper:

> A. Gachahi et al., "Revolutionizing Wireless Communication with Single Layer Capacitor-Based Antenna Technology", arXiv:2403.01170.

---

## 🌐 Connect with Us

For latest updates, development kits, and more advanced simulations:

* Visit [**EMF Systems Official Website**](https://emfsystems.africa)
* Follow us on [LinkedIn](https://www.linkedin.com/emf systems)

---

## 📷 Example Outputs

*Example radiation pattern from a 4x4 SLC compound array:*
