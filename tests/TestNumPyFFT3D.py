import sys
import numpy as np

#-------------- initialize ------------
print("Initializing")
np.set_printoptions(precision=10)
data_init = [
0.183471406e+0,0.623968915e+0,0.731257661e+0,0.997228140e+0,0.961913696e+0,
0.792673860e-1,0.429684069e+0,0.290531312e+0,0.453270921e+0,0.199228629e+0,
0.754931905e-1,0.226924328e+0,0.936407886e+0,0.979392715e+0,0.464957186e+0,
0.742653949e+0,0.368019859e+0,0.885231224e+0,0.406191773e+0,0.653096157e+0,
0.567929080e-1,0.568028857e+0,0.144986181e+0,0.466158777e+0,0.573327733e+0,
0.136324723e+0,0.819010407e+0,0.271218167e+0,0.626224101e+0,0.398109186e-1,
0.860031651e+0,0.338153865e+0,0.688078522e+0,0.564682952e+0,0.222924187e+0,
0.306816449e+0,0.316316038e+0,0.640568415e+0,0.702342408e+0,0.632135481e+0,
0.649402777e+0,0.647100865e+0,0.370402133e+0,0.691313864e+0,0.447870566e+0,
0.757298851e+0,0.586173682e+0,0.766745717e-1,0.504185402e+0,0.812016428e+0,
0.217988206e+0,0.273487202e+0,0.937672578e+0,0.570540523e+0,0.409071185e+0,
0.391548274e-1,0.663478965e+0,0.260755447e+0,0.503943226e+0,0.979481790e+0]
data_k_answer = [
 30.0601362322000+0.000000000000000e+0j    , 0.353642310400000-0.656637882635999j,
 1.84441281060000-2.74233574840000j        ,-2.46150775700102+4.133457522749440e-2j,
 0.817180262600000+0.000000000000000e+0j   ,-0.825032729800000+0.726831889225593j,
 1.84441281060000+2.74233574840000j        , 0.732077908401022-1.00081656417251j,
-0.458084677221565+0.153852624198044j      ,-0.144717578762292-2.28803500866730j,
 0.705009791908085-2.01933903816772j       , 2.964764771912920e-2-2.25308567091476j,
-0.973844201363205-2.06176730216391j       ,-0.997469264485567-0.602564694796768j,
 1.40462480063246+1.01369298185703j        ,1.75359066116106+0.326376682455054j,
-1.44138430512843-1.31255361296267j        ,-0.744097894016966-0.480832819344072j,
-2.15508603929787-1.93166708304686j        ,-0.438393767508181-1.10483993164215j,
 0.995576356313205+0.780345054198534j      ,-3.86510390713435+4.55728094860600j,
 2.709703615732217e-2-2.321947897161802e-2j,0.208459929713295+2.48927791839840j,
-1.44138430512843+1.31255361296267j        ,1.68993207144876-0.218131499149520j,
 2.709703615732217e-2+2.321947897161802e-2j,-2.04216524355430+2.34581169104130j,
 0.995576356313205-0.780345054198534j      ,2.469425255987412e-2+1.01595227927158j,
-2.15508603929787+1.93166708304686j        ,-0.744464511100393-5.438522180051475e-2j,
-0.458084677221565-0.153852624198044j      ,-0.713266212819504+1.64663971056941j,
 1.40462480063246-1.01369298185703j        ,-2.32489174723348-1.41241859006073j,
-0.973844201363205+2.06176730216391j       ,0.857829657610043-1.36198877135770j,
 0.705009791908085+2.01933903816772j       ,-0.231601465597141+0.142526551270715j]

data_init = np.reshape(data_init, (5,4,3))
data_k_answer = np.reshape(data_k_answer, (5,4,2))
#---------------- Forward --------------------
print("Running FFT 3D")
print("If error is less than 1.0e-7, it is ok!")
data_k = np.fft.rfftn(data_init)
error = np.max(np.absolute(data_k - data_k_answer))
print("FFT Forward Error: ", error)
if(np.isnan(error) or error > 1e-7):
    sys.exit(-1)

#--------------- Backward --------------------
data_r = np.fft.irfftn(data_k_answer, (5,4,3))
error = np.max(np.absolute(data_r - data_init))
print("FFT Backward Error: ", error)
if(np.isnan(error) or error > 1e-7):
    sys.exit(-1)

"""
#--------------- Test with large array --------------------
data_init = np.random.uniform(-1.0,1.0, size=(89,101,119))
data_k = np.fft.rfftn(data_init)
data_r = np.fft.irfftn(data_k, (data_init.shape))
error = np.max(np.absolute(data_r - data_init))
print("Test with lage array, Error: ", error)
if(np.isnan(error) or error > 1e-7):
    sys.exit(-1)
"""
