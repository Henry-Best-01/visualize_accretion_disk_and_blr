A simple Streamlit app showing how Amoeba can produce flux distributions of the active galactic nuclei (AGN) accretion disk and broad line region (BLR). Various sliders control the typical accretion disk parameters and our viewing conditions. Applying the "GR" toggle will apply relativistic corrections via pre-calculated raytraces computed using Sim5. 

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://gdp-dashboard-template.streamlit.app/)

### How to run it on your own machine

1. Install the requirements

   ```
   $ pip install -r requirements.txt
   ```

2. Run the app

   ```
   $ streamlit run visualize_AD_and_BLR.py
   ```
