# 🔬 CT Reconstruction Lab

A web-based interactive tool for simulating Computed Tomography (CT) scans and reconstruction algorithms.

## 🌟 Overview
This project provides an interactive simulation of the CT imaging pipeline. It allows users to select and configure digital phantoms, simulate X-ray projections using the Radon Transform, model realistic Poisson noise (Beer-Lambert Law), and compare various analytic and iterative reconstruction algorithms.

## 🚀 Live Demo
You can try the app directly in your browser here: https://ct-reconstruction-lab-orcqs9e2rztuqn8eheyv6p.streamlit.app/

## 🛠️ Features
- **Phantom Generation:** Shepp-Logan, Point Source, and custom geometric shapes.
- **Sinogram Simulation:** Adjustable projection angles and X-ray intensity (Dose).
- **Noise Modeling:** Realistic Poisson noise simulation.
- **Reconstruction Algorithms:**
    - Analytic: Back Projection (BP) and Filtered Back Projection (FBP).
    - Iterative: SIRT and SART.

## 💻 Local Installation
If you want to run this locally:
1. Clone the repo:
```bash
git clone https://github.com/sadragh23/ct-reconstruction-lab.git
```
2. Install dependencies: 
```bash
pip install -r requirements.txt
```
3. Run the app: 
```bash 
streamlit run streamlit_app.py
```