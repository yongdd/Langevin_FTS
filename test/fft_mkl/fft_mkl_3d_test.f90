!--------------------------------------------------------------------
! The main program begins here
program block_copolymer_3d

  use fft
  implicit none

  integer, parameter :: II = 3
  integer, parameter :: JJ = 4
  integer, parameter :: KK = 5

  integer, parameter :: BOX_DIM = 3

  real(kind=8) :: data_init(1:II,1:JJ,1:KK)
  real(kind=8) :: data_r(1:II,1:JJ,1:KK)
  complex(kind=8) :: data_k(1:II/2+1,1:JJ,1:KK)
  complex(kind=8) :: data_k_answer(1:II/2+1,1:JJ,1:KK)
!-------------- initialize ------------
  write(*,*) "Initializing."

  call fft_initialize(BOX_DIM, &
    II, JJ, KK)
!---------------- run --------------------
  write(*,*) "Running FFT MKL 3D."

  data_init = reshape((/0.183471406D+00,0.623968915D+00,0.731257661D+00,0.997228140D+00,0.961913696D+00, &
    0.792673860D-01,0.429684069D+00,0.290531312D+00,0.453270921D+00,0.199228629D+00, &
    0.754931905D-01,0.226924328D+00,0.936407886D+00,0.979392715D+00,0.464957186D+00, &
    0.742653949D+00,0.368019859D+00,0.885231224D+00,0.406191773D+00,0.653096157D+00, &
    0.567929080D-01,0.568028857D+00,0.144986181D+00,0.466158777D+00,0.573327733D+00, &
    0.136324723D+00,0.819010407D+00,0.271218167D+00,0.626224101D+00,0.398109186D-01, &
    0.860031651D+00,0.338153865D+00,0.688078522D+00,0.564682952D+00,0.222924187D+00, &
    0.306816449D+00,0.316316038D+00,0.640568415D+00,0.702342408D+00,0.632135481D+00, &
    0.649402777D+00,0.647100865D+00,0.370402133D+00,0.691313864D+00,0.447870566D+00, &
    0.757298851D+00,0.586173682D+00,0.766745717D-01,0.504185402D+00,0.812016428D+00, &
    0.217988206D+00,0.273487202D+00,0.937672578D+00,0.570540523D+00,0.409071185D+00, &
    0.391548274D-01,0.663478965D+00,0.260755447D+00,0.503943226D+00,0.979481790D+00/),(/II,JJ,KK/))

  data_k_answer = reshape((/(30.0601362322000,0.000000000000000E+000), (0.353642310400000,-0.656637882635999), &
    (1.84441281060000,-2.74233574840000), (-2.46150775700102,4.133457522749440E-002), &
    (0.817180262600000,0.000000000000000E+000), (-0.825032729800000,0.726831889225593), &
    (1.84441281060000,2.74233574840000), (0.732077908401022,-1.00081656417251), &
    (-0.458084677221565,0.153852624198044), (-0.144717578762292,-2.28803500866730), &
    (0.705009791908085,-2.01933903816772), (2.964764771912920E-002,-2.25308567091476), &
    (-0.973844201363205,-2.06176730216391), (-0.997469264485567,-0.602564694796768), &
    (1.40462480063246,1.01369298185703), (1.75359066116106,0.326376682455054), &
    (-1.44138430512843,-1.31255361296267), (-0.744097894016966,-0.480832819344072), &
    (-2.15508603929787,-1.93166708304686), (-0.438393767508181,-1.10483993164215), &
    (0.995576356313205,0.780345054198534), (-3.86510390713435,4.55728094860600), &
    (2.709703615732217E-002,-2.321947897161802E-002), (0.208459929713295,2.48927791839840), &
    (-1.44138430512843,1.31255361296267), (1.68993207144876,-0.218131499149520), &
    (2.709703615732217E-002,2.321947897161802E-002), (-2.04216524355430,2.34581169104130), &
    (0.995576356313205,-0.780345054198534), (2.469425255987412E-002,1.01595227927158), &
    (-2.15508603929787,1.93166708304686), (-0.744464511100393,-5.438522180051475E-002), &
    (-0.458084677221565,-0.153852624198044), (-0.713266212819504,1.64663971056941), &
    (1.40462480063246,-1.01369298185703), (-2.32489174723348,-1.41241859006073), &
    (-0.973844201363205,2.06176730216391), (0.857829657610043,-1.36198877135770), &
    (0.705009791908085,2.01933903816772), (-0.231601465597141,0.142526551270715)/),(/II/2+1,JJ,KK/))

  call fft_forward(data_init,data_k)
  call fft_backward(data_k_answer,data_r)
  
!--------------- check --------------------
  write(*,*) "Checking."
  write(*,*) "If error is less than 1.0e-6, it is ok!"

  write(*,*) "fft_forward error:",  sqrt(maxval(abs(data_k-data_k_answer)**2))
  write(*,*) "fft_backward error:", sqrt(maxval(abs(data_r-data_init)**2))

!------------- finalize -------------
  write(*,*) "Finalize."
! finalize pseudo_spectral module
  call fft_finalize()

end program block_copolymer_3d
