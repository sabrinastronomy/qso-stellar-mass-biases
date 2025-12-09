import numpy as np
import matplotlib.pyplot as plt


def L_to_M(L,z):   # from erg/s to Muv (Fontanot+12, see Ni+19)
   #Values taken from Gnedin's cosmological calculator, assuming Omega0=0.3278, H0=69.7km/s/Mpc
   if z==6.5:
       dlMpc=64863.7 #Luminosity distance, CHECK THESE
       dm=49.06 #Distance modulus, CHECK THESE

   if z==7:
       dlMpc=70878.29 #Luminosity distance
       dm=49.25 #Distance modulus
       dmK=46.99 #Distance modulus + k-correction
   if z==8:
       dlMpc=82688.79
       dm=49.59
       dmK=47.2
   if z==9:
       dlMpc=94637.09
       dm=49.88
       dmK=47.38
   if z==10:
       dlMpc=106715.22
       dm=50.14
       dmK=47.54
   if z==11:
       dlMpc=118912.91
       dm=50.38
       dmK=47.68
   if z==12:
       dlMpc=131198
       dm=50.59
       dmK=dm-2.78

   dl=(dlMpc*1e6*3.086e+18) #cm
   M = -2.5 * np.log10(L/(4*np.pi*dl**2))- 48.60-dm
   #No k-correction required for LUV -> MUV
   return M

z = 6.5

f356_mags = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/full_6_JWST.NIRCAM.F356W_mag_all_order_luminous_BH.npy").flatten()
abs_lum = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/full_6_FAKE.FAKE.1500_mag_all_order_luminous_BH.npy").flatten()
maddie_uv = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/maddie_lum_files/lum_dust_6p5_1500.npy")

abs_mags = L_to_M(abs_lum, z)
abs_mags_maddie = L_to_M(maddie_uv, z)

plt.scatter(f356_mags, abs_mags)
plt.xlabel("f356")
plt.ylabel("absolute uv")
plt.savefig("../image_pngs/comp.png")
plt.close()

plt.plot(abs_mags,  linestyle='None', marker='o')
plt.ylabel("UV mag")
plt.xlabel("galaxy number")
plt.savefig("../image_pngs/uv_mag_mine.png")

plt.close()

plt.plot(abs_mags_maddie,  linestyle='None', marker='o')
plt.ylabel("UV mag (maddie)")
plt.xlabel("galaxy number")
plt.savefig("../image_pngs/uv_mag_maddie.png")

plt.close()

print(len(f356_mags))