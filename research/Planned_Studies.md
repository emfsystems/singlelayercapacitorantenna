# 📅 Planned Studies and Simulation Extensions

This document outlines the **future simulation studies** and **planned expansions** of the Single Layer Capacitor-Based Antenna project.

These studies will build upon the core capabilities of the `SingleLayerCapacitorAntenna_new` class to demonstrate real-world integration across multiple wireless technologies and scientific applications.

**Note:** This list is not exhaustive. As the project evolves, additional areas of exploration may be added or existing areas may expand based on new opportunities, industry trends, or community contributions.

---

## 🛶 1. WiFi Network Integration
- Objective: Model SLC antennas as access point antennas in WiFi 5, WiFi 6, and WiFi 7 networks, including **Mesh WiFi Networks** designed for **long-range wireless communication** that can rival GSM cellular networks.
- Goals:
  - Analyze indoor and outdoor propagation characteristics.
  - Evaluate mesh node performance, range extension, and redundancy.
  - Compare against conventional PCB trace antennas.
- Outputs:
  - Radiation patterns
  - Coverage maps
  - Throughput simulations (using MATLAB WLAN Toolbox)


## 🌌 2. Satellite Communication Links
- Objective: Analyze SLC antenna performance in satellite uplink and downlink scenarios (LEO/MEO/GEO) and extend simulations to **Deep Space Communication** systems (e.g., Earth-Moon, Earth-Mars links).
- Goals:
  - Simulate link budgets and antenna gain requirements for Earth-orbit and deep space missions.
  - Study beam shaping, polarization control, and delay compensation for long-range links.
- Outputs:
  - Link budget calculation scripts
  - Radiation pattern shaping examples
  - Deep space communication feasibility studies


## 🌐 3. 5G Base Station Arrays
- Objective: Simulate SLC antennas in 5G massive MIMO applications.
- Goals:
  - Phased array configurations
  - Beamforming studies at 3.5 GHz and mmWave bands (28 GHz, 39 GHz)
  - Study multi-user MIMO and dynamic beam steering.
- Outputs:
  - Beamforming patterns
  - S-parameter matrices for subarrays


## 💻 4. IoT Gateway Antennas
- Objective: Model antennas for Sub-GHz and 2.4 GHz IoT devices and gateways.
- Goals:
  - Analyze long-range coverage in dense and rural environments.
  - Optimize energy efficiency in battery-powered or energy-harvesting devices.
- Outputs:
  - Coverage simulation maps
  - Efficiency vs. distance charts


## 🌊 5. LPWAN (LoRa, NB-IoT) Integration
- Objective: Test SLC antennas for LPWAN applications.
- Goals:
  - LoRa gateway efficiency and long-range node performance.
  - Low-power optimizations for NB-IoT devices.
- Outputs:
  - Range simulation studies
  - Packet delivery ratio predictions


## 📷 6. Imaging, Radar Sensors, and Radio Astronomy
- Objective: Investigate the use of SLC antennas in scientific sensing, imaging, and astronomical observation systems.
- Goals:
  - Simulate active and passive radar systems using SLC antenna arrays.
  - Explore imaging techniques such as synthetic aperture radar (SAR).
  - Model radio telescope arrays for deep-space signal reception.
- Outputs:
  - Sensor system simulation results
  - Imaging resolution and sensitivity studies
  - Radio astronomy application demonstrations


---

# 🔄 Update Process
- Planned studies will evolve as new simulation capabilities are added.
- Each study will eventually include:
  - A new example script under `/examples/`
  - Results published under `/research/datasets/` and `/research/papers/`

This living document will grow as the project expands.

---

# 🔹 End of Planned Studies
