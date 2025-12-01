
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

# Carica tutti i dati
data = np.load('atmospheric_data.npz')


# Crea il grafico
alt_array=np.logspace(0, 6, 100)  # Altitudine da 0 a 1000 km
msise90_low_interp = interp.interp1d(data['msise90_low_altitude'], np.log(data['msise90_low_density']), kind='linear', fill_value="extrapolate")
msise90_mean_interp = interp.interp1d(data['msise90_mean_altitude'], np.log(data['msise90_mean_density']), kind='linear', fill_value="extrapolate")
msise90_high_interp = interp.interp1d(data['msise90_high_altitude'], np.log(data['msise90_high_density']), kind='linear', fill_value="extrapolate")
plt.plot(np.exp(msise90_low_interp(alt_array)),  alt_array*1e-3, label='MSISE90 Low', color='blue',marker='.',ls='None')
plt.plot(np.exp(msise90_mean_interp(alt_array)), alt_array*1e-3, label='MSISE90 Mean', color='orange',marker='.',ls='None')
plt.plot(np.exp(msise90_high_interp(alt_array)), alt_array*1e-3, label='MSISE90 High', color='green',marker='.',ls='None')
plt.xscale('log')
plt.ylabel('Altitudine (km)')   
plt.xlabel('Densità (kg/m³)')   
plt.title('Densità Atmosferica vs Altitudine (MSISE90)')
plt.legend()
plt.grid(True)
plt.show()