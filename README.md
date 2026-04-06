# CT Reconstruction Lab

Interactive CT reconstruction lab for exploring Radon transform, sinograms, noise effects, and iterative reconstruction algorithms (BP, FBP, SART, SIRT) using Python and Streamlit.

## Features
- Multiple phantom types (Shepp-Logan, Square, Point Source, Two Circles)
- Sinogram generation using Radon transform
- Optional Poisson noise simulation
- Reconstruction methods:
  - Back Projection (BP)
  - Filtered Back Projection (FBP)
  - SART
  - SIRT

## Tech Stack
- Python
- NumPy
- scikit-image
- Matplotlib
- Streamlit
## MATLAB
A MATLAB implementation of the classical ART algorithm is provided in the `/ART MATLAB` folder.

## How to Run
```bash
pip install -r requirements.txt
streamlit run main.py
```
## Live Demo
