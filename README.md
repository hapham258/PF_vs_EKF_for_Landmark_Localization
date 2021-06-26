# PF_vs_EKF_for_Landmark_Localization
Performance Comparison of Particle Filter and Extended Kalman Filter for Landmark-based Localization, implemented in MATLAB.

## Summary
This work is my assignment in the course Artificial Intelligence, taught by Dr. Viet-Cuong Pham.

## Results
<p align="center">
  <img src="ekf.svg" width="400" alt="accessibility text">
  <img src="pf_best.svg" width="400" alt="accessibility text">
</p>
<p align="center">
  <img src="pf_avg.svg" width="400" alt="accessibility text">
  <img src="pf_worst.svg" width="400" alt="accessibility text">
</p>

-- | x (m) | y (m) | theta (rad)
-- | -- | -- | --
error_odo | 0.612770 | 0.815380 | 0.434710
error_max | 0.076751 | 0.086025 | 0.012091
error_avg | 0.077178 | 0.106180 | 0.012687
error_min | 0.090608 | 0.130630 | 0.015909
error_ekf | 0.136460 | 0.108730 | 0.302300
