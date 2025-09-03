Extreme Precipitation Return Levels (BHM × GHCNd) 🌧️

Generate high-resolution maps of extreme-precipitation return levels using a Bayesian Hierarchical Model (BHM) trained on the Global Historical Climatology Network – daily (GHCNd) dataset.

What this repo does

Modeling:

Location (μ) and scale (σ) are modeled via Gaussian Processes.

Shape (ξ) is treated as spatially homogeneous due to its high estimation uncertainty.

Covariates: Longitude, Latitude, Elevation.

A model-selection step prunes unnecessary covariates (described in the upcoming paper).

Outputs: Annual (and optionally seasonal) T-year return level rasters and figures.

✅ Annual extreme-precipitation dataset is included.
📬 Seasonal datasets & additional return periods: contact atakallou@crimson.ua.edu




<img width="1275" height="610" alt="image" src="https://github.com/user-attachments/assets/c2307b43-ee34-4c11-91b2-6c928dd9e178" />
<img width="1275" height="614" alt="image" src="https://github.com/user-attachments/assets/237673c5-70ff-4a65-af39-f9fc6f662eef" />
