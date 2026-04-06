import streamlit as st
import numpy as np
from skimage.data import shepp_logan_phantom
from skimage.draw import disk
from skimage.transform import resize,radon, iradon, iradon_sart
import matplotlib.pyplot as plt

# Set layout to wide to accommodate three columns comfortably
st.set_page_config(layout="wide", page_title="CT Reconstruction Lab")

st.title("🔬 Computed Tomography Reconstruction")
st.markdown("---")

# Create the three main columns
col1, col2, col3 = st.columns(3)

# --- COLUMN 1: ORIGINAL PHANTOM ---
with col1:
    st.subheader("Original Image")
    
    with st.expander("Phantom Settings", expanded=True):
        phantom_type = st.selectbox("Select Phantom Type", 
                                    ["Shepp-Logan", "Simple Square", "Point Source", "Two Circles"])
        N = st.select_slider("Image Resolution", options=[64,128,256], value=128)

        # create the phantoms
        if phantom_type == "Shepp-Logan":
            img_true = resize(shepp_logan_phantom(), (N, N))
        elif phantom_type == "Simple Square":
            img_true = np.zeros((N, N))
            start = N // 3
            end = start + (N // 6) 
            img_true[start:end, start:end] = 1.0
        elif phantom_type == "Point Source":
            img_true = np.zeros((N, N))
            img_true[N // 2 - 20, N // 2 - 25] = 1.0
        elif phantom_type == "Two Circles":
            img_true = np.zeros((N, N))
            # Circle 1: Slightly top-left of center
            rr1, cc1 = disk((N // 2 - N // 8, N // 2 - N // 8), N // 12)
            img_true[rr1, cc1] = 1.0
            # Circle 2: Slightly bottom-right of center
            rr2, cc2 = disk((N // 2 + N // 8, N // 2 + N // 8), N // 15)
            img_true[rr2, cc2] = 0.7  # Different intensity to see the contrast
        
    # plot the image
    st.subheader(f"Generated {phantom_type} Phantom")
    fig, ax = plt.subplots()
    ax.imshow(img_true, cmap='gray')
    ax.axis('off')
    st.pyplot(fig)

# --- COLUMN 2: SINOGRAM (RADON TRANSFORM) ---
with col2:
    st.subheader("Sinogram")
    
    with st.expander("Projections Settings", expanded=True):
        angles = st.number_input("Number of Angles", min_value=1, max_value=360, value=180)
        step = st.number_input("Step Size (Degrees)", min_value=0.1, max_value=10.0, value=1.0)
        
        # Add noise
        add_noise = st.checkbox("Add Poisson Noise")
        # angles
        theta = np.arange(0,angles,step)
        # create sinogram
        sinogram = radon(img_true, theta=theta)

        if add_noise:
            I0 = st.select_slider("X-ray Intensity (Dose)", [10, 100, 1000, 10000, 100000])
            # 1. Normalize
            m = np.max(sinogram)
            if m > 0:
                P = sinogram / m
            else:
                P = sinogram
            
            # 2. Physics: Convert to photon counts (Beer-Lambert Law)
            counts = I0 * np.exp(-P)
            
            # 3. Add Poisson noise (Random photon arrival)
            noisy_counts = np.random.poisson(counts).astype(float)
            
            # 4. Clean up zeros (you can't take log of 0)
            noisy_counts[noisy_counts <= 0] = 1e-6
            
            # 5. Convert back to sinogram (Linear Attenuation)
            sinogram = np.log(I0 / noisy_counts)
            sinogram = sinogram * m
        #save to session state
        st.session_state.sinogram = sinogram
        st.session_state.theta = theta

    fig_sino, ax_sino = plt.subplots(figsize=(6, 4))
    ax_sino.imshow(
        sinogram, 
        cmap='gray', 
        interpolation='nearest',
        aspect='auto', #Let the plot stretch so it's not a tiny square
        extent=(theta.min(), theta.max(), 0, sinogram.shape[0])
    )
    st.pyplot(fig_sino)

# --- COLUMN 3: RECONSTRUCTION RESULT ---
with col3:
    st.subheader("Reconstruction")

    sinogram = st.session_state.sinogram
    theta = st.session_state.theta
    N = sinogram.shape[0]

    with st.expander("Algorithm Settings", expanded=True):
        algo = st.selectbox("Method",
                            ["Back Projection(BP)", "Filtered Back Projection(FBP)", "SART", "SIRT"])

        if algo == "Filtered Back Projection(FBP)":
            filter_type = st.select_slider("Filter", options=["ramp", "shepp-logan", "cosine"])
        if algo == "SIRT":
            n_iter = st.number_input("Iterations(1-200)", min_value=1, max_value=200, value=100)
            relaxation = st.select_slider("Relaxation",
                                        [round(x,2) for x in np.arange(0.05, 1.05, 0.05)],
                                        value=0.2)
        if algo == "SART":
            n_iter = st.number_input("Iterations(1-200)", min_value=1, max_value=30, value=10)
            relaxation = st.select_slider("Relaxation",
                                        [round(x,2) for x in np.arange(0.01, 1.0, 0.01)],
                                        value=0.2)



        # Execute BP
        if algo == "Back Projection(BP)":
            rec = iradon(sinogram, theta=theta, filter_name=None)
        # Execute FBP
        if algo == "Filtered Back Projection(FBP)":
            rec = iradon(sinogram, theta=theta, filter_name=filter_type)
        # Ececute SIRT (Simultaneous Iterative Reconstruction Technique)
        if algo == "SIRT":
            # compares the original sinogram with estimated sinogram and update the image
            
            #start with zero image
            rec = np.zeros((N,N))
            n_projections = len(theta)

            progress_bar = st.progress(0) #progress bar
            for i in range(n_iter):
                # forward project of our current guess
                sino_guess = radon(rec, theta=theta, circle=True)

                #find the Error (difference between measured and guess)
                error_sino = sinogram - sino_guess

                # Back project the error
                error_rec = iradon(error_sino, theta=theta, filter_name=None, circle=True)

                # update image
                rec += (relaxation / n_projections) * error_rec

                #progress bar
                progress_bar.progress((i + 1) / n_iter)

        # Execute SART    
        if algo == "SART":
            rec = np.zeros((N,N))
            for j in range(n_iter):
                #SART updates the image per-angle.
                rec = iradon_sart(sinogram, theta=theta, image=rec, relaxation=relaxation)


    # Plot the Results
    st.subheader("reconstructed image")
    fig_rec, ax_rec = plt.subplots()
    ax_rec.imshow(rec, cmap='gray')
    ax_rec.axis('off')
    st.pyplot(fig_rec)
    
    
    # Information (Metrics)

    #Normalize the images
    img_ref = (img_true - img_true.min()) / (img_true.max() - img_true.min() + 1e-9)
    img_res = (rec - rec.min()) / (rec.max() - rec.min() + 1e-9)
    
    #calculate mean squared error
    mse_value = np.mean((img_ref - img_res) ** 2)

    st.info("Reconstruction Metrics")
    st.write(f"**MSE:** {mse_value:.4f}")
    
    if st.error: # Example of where you'd show errors
        pass 


# Extra row for analysis.

#st.markdown("---")

#st.subheader("Analysis & Residuals")
#bot_col1, bot_col2 = st.columns(2)

#with bot_col1:
#    st.write("Difference Map (Original vs Recon)")
#    st.image("https://via.placeholder.com/400x200", caption="Residuals")

#with bot_col2:
#    st.write("Profile Comparison")
#    # You can put a line chart here later
#    st.line_chart([0, 1, 2, 3, 2, 1, 0])
