import os
import sys
import copy

_sacrum = ("sacrum")
_pelvis = ("pelvis")
#
_pelvis_l = ("pelvis_l")
_pelvis_r = ("pelvis_r")
#
_femur = ("femur")
_femur_l = ("femur_l")
_femur_r = ("femur_r")
#
_adductor=("adductor")
_adductor_l=("adductor_l")
_adductor_r=("adductor_r")
#
_biceps_femoris=("biceps_femoris")
_biceps_femoris_l=("biceps_femoris_l")
_biceps_femoris_r=("biceps_femoris_r")
#
_gluteus_maximus=("gluteus_maximus")
_gluteus_maximus_l=("gluteus_maximus_l")
_gluteus_maximus_r=("gluteus_maximus_r")
#
_gluteus_medius=("gluteus_medius")
_gluteus_medius_l=("gluteus_medius_l")
_gluteus_medius_r=("gluteus_medius_r")
#
_gluteus_minimus=("gluteus_minimus")
_gluteus_minimus_l=("gluteus_minimus_l")
_gluteus_minimus_r=("gluteus_minimus_r")
#
_gracilis=("gracilis")
_gracilis_l=("gracilis_l")
_gracilis_r=("gracilis_r")
#
_iliacus=("iliacus")
_iliacus_l=("iliacus_l")
_iliacus_r=("iliacus_r")
#
_obturator_externus=("obturator_externus")
_obturator_externus_l=("obturator_externus_l")
_obturator_externus_r=("obturator_externus_r")
#
_obturator_internus=("obturator_internus")
_obturator_internus_l=("obturator_internus_l")
_obturator_internus_r=("obturator_internus_r")
#
_pectineus=("pectineus")
_pectineus_l=("pectineus_l")
_pectineus_r=("pectineus_r")
#
_piriformis=("piriformis")
_piriformis_l=("piriformis_l")
_piriformis_r=("piriformis_r")
#
_psoas_major=("psoas_major")
_psoas_major_l=("psoas_major_l")
_psoas_major_r=("psoas_major_r")
#
_rectus_femoris=("rectus_femoris")
_rectus_femoris_l=("rectus_femoris_l")
_rectus_femoris_r=("rectus_femoris_r")
#
_sartorius=("sartorius")
_sartorius_l=("sartorius_l")
_sartorius_r=("sartorius_r")
#
_semimembranosus=("semimembranosus")
_semimembranosus_l=("semimembranosus_l")
_semimembranosus_r=("semimembranosus_r")
#
_semitendinosus=("semitendinosus")
_semitendinosus_l=("semitendinosus_l")
_semitendinosus_r=("semitendinosus_r")
#
_tensor_fasciae_latae=("tensor_fasciae_latae")
_tensor_fasciae_latae_l=("tensor_fasciae_latae_l")
_tensor_fasciae_latae_r=("tensor_fasciae_latae_r")
#
_vastus_lateralis_intermedius=("vastus_lateralis_intermedius")
_vastus_lateralis_intermedius_l=("vastus_lateralis_intermedius_l")
_vastus_lateralis_intermedius_r=("vastus_lateralis_intermedius_r")
#
_vastus_medialis=("vastus_medialis")
_vastus_medialis_l=("vastus_medialis_l")
_vastus_medialis_r=("vastus_medialis_r")
#
_tibia=("tibia")
_tibia_l=("tibia_l")
_tibia_r=("tibia_r")
#
_fibula=("fibula")
_fibula_l=("fibula_l")
_fibula_r=("fibula_r")
#
_anterior_comp=("anterior_comp")
_anterior_comp_l=("anterior_comp_l")
_anterior_comp_r=("anterior_comp_r")
#
_lateral_comp=("lateral_comp")
_lateral_comp_l=("lateral_comp_l")
_lateral_comp_r=("lateral_comp_r")
#
_deep_posterior_comp=("deep_posterior_comp")
_deep_posterior_comp_l=("deep_posterior_comp_l")
_deep_posterior_comp_r=("deep_posterior_comp_r")
#
_superficial_posterior_comp=("superficial_posterior_comp")
_superficial_posterior_comp_l=("superficial_posterior_comp_l")
_superficial_posterior_comp_r=("superficial_posterior_comp_r")
#
#
_quadratus_femoris_muscle=("quadratus_femoris_muscle")
_quadratus_femoris_muscle_l=("quadratus_femoris_muscle_l")
_quadratus_femoris_muscle_r=("quadratus_femoris_muscle_r")
#
_vastus_intermedius_muscle=("vastus_intermedius_muscle")
_vastus_intermedius_muscle_l=("vastus_intermedius_muscle_l")
_vastus_intermedius_muscle_r=("vastus_intermedius_muscle_r")
#
#
_talus=("talus")
_talus_l=("talus_l")
_talus_r=("talus_r")
#
_talus_center=("talus_center")
_talus_center_l=("lt_talus_center")
_talus_center_r=("rt_talus_center")
#
_calcaneal=("calcaneal")
_calcaneal_l=("calcaneal_l")
_calcaneal_r=("calcaneal_r")
#
_calcaneal_center=("calcaneal_center")
_calcaneal_center_l=("lt_calcaneal_center")
_calcaneal_center_r=("rt_calcaneal_center")
#
_navicular=("navicular")
_navicular_l=("navicular_l")
_navicular_r=("navicular_r")
#
_navicular_center=("navicular_center")
_navicular_center_l=("lt_navicular_center")
_navicular_center_r=("rt_navicular_center")
#
_medial_cuneiform=("medial_cuneiform")
_medial_cuneiform_l=("medial_cuneiform_l")
_medial_cuneiform_r=("medial_cuneiform_r")
#
_medial_cuneiform_center=("medial_cuneiform_center")
_medial_cuneiform_center_l=("lt_medial_cuneiform_center")
_medial_cuneiform_center_r=("rt_medial_cuneiform_center")
#
_intermediate_cuneiform=("intermediate_cuneiform")
_intermediate_cuneiform_l=("intermediate_cuneiform_l")
_intermediate_cuneiform_r=("intermediate_cuneiform_r")
#
_intermediate_cuneiform_center=("intermediate_cuneiform_center")
_intermediate_cuneiform_center_l=("lt_intermediate_cuneiform_center")
_intermediate_cuneiform_center_r=("rt_intermediate_cuneiform_center")
#
_lateral_cuneiform=("lateral_cuneiform")
_lateral_cuneiform_l=("lateral_cuneiform_l")
_lateral_cuneiform_r=("lateral_cuneiform_r")
#
_lateral_cuneiform_center=("lateral_cuneiform_center")
_lateral_cuneiform_center_l=("lt_lateral_cuneiform_center")
_lateral_cuneiform_center_r=("rt_lateral_cuneiform_center")
#
_cuboid=("cuboid")
_cuboid_l=("cuboid_l")
_cuboid_r=("cuboid_r")
#
_cuboid_center=("cuboid_center")
_cuboid_center_l=("lt_cuboid_center")
_cuboid_center_r=("rt_cuboid_center")
#
_foot_bone=("foot_bone")
_foot_bone_l=("foot_bone_l")
_foot_bone_r=("foot_bone_r")
#
_patella=("patella")
_patella_l=("patella_l")
_patella_r=("patella_r")
#
_skin=("skin")
#
_origin=("origin")
_origin_point=("origin_point")
#
_skin_leg=("skin_leg")
_skin_upperleg=("skin_upperleg")
_skin_lowerleg=("skin_lowerleg")
#
_skin_pelvis=("skin_pelvis")
_skin_pelvis_l=("skin_pelvis_l")
_skin_pelvis_r=("skin_pelvis_r")
_skin_femur=("skin_femur")
_skin_femur_l=("skin_femur_l")
_skin_femur_r=("skin_femur_r")
_skin_tibia=("skin_tibia")
_skin_tibia_l=("skin_tibia_l")
_skin_tibia_r=("skin_tibia_r")
_skin_foot_bone=("skin_foot_bone")
_skin_foot_bone_l=("skin_foot_bone_l")
_skin_foot_bone_r=("skin_foot_bone_r")
#
_pubic_tubercle = ("pubic_tubercle")
#
_anterior_sacral_slope = ("anterior_sacral_slope")
_posterior_sacral_slope = ("posterior_sacral_slope")
_ASIS = ("ASIS")
_rt_ASIS = ("rt_ASIS")
_lt_ASIS = ("lt_ASIS")
_hip_center = ("hip_center")
_rt_hip_center = ("rt_hip_center")
_lt_hip_center = ("lt_hip_center")
#
_knee_center=("knee_center")
_rt_knee_center=("rt_knee_center")
_lt_knee_center=("lt_knee_center")
_ankle_center=("ankle_center")
_rt_ankle_center=("rt_ankle_center")
_lt_ankle_center=("lt_ankle_center")
_prox_femur=("prox_femur")
_rt_prox_femur=("rt_prox_femur")
_lt_prox_femur=("lt_prox_femur")
_late_retrocondy=("late_retrocondy")
_rt_late_retrocondy=("rt_late_retrocondy")
_med_retrocondy=("med_retrocondy")
_rt_med_retrocondy=("rt_med_retrocondy")
_lt_late_retrocondy=("lt_late_retrocondy")
_lt_med_retrocondy=("lt_med_retrocondy")
_late_retrotibia=("late_retrotibia")
_rt_late_retrotibia=("rt_late_retrotibia")
_med_retrotibia=("med_retrotibia")
_rt_med_retrotibia=("rt_med_retrotibia")
_lt_late_retrotibia=("lt_late_retrotibia")
_lt_med_retrotibia=("lt_med_retrotibia")
_late_malleolus=("late_malleolus")
_rt_late_malleolus=("rt_late_malleolus")
_med_malleolus=("med_malleolus")
_rt_med_malleolus=("rt_med_malleolus")
_lt_late_malleolus=("lt_late_malleolus")
_lt_med_malleolus=("lt_med_malleolus")
#
_toe = ("toe")
_toe_l=("lt_toe")
_toe_r=("rt_toe")
#
_name2id=OrderedDict({})
#
_name2id[_pelvis_l] = 1
_name2id[_pelvis_r] = 2
_name2id[_femur_l] = 3
_name2id[_femur_r] = 4
_name2id[_adductor_l] = 5
_name2id[_adductor_r] = 6
_name2id[_biceps_femoris_l] = 7
_name2id[_biceps_femoris_r] = 8
_name2id[_gluteus_maximus_l] = 9
_name2id[_gluteus_maximus_r] = 10
_name2id[_gluteus_medius_l] = 11
_name2id[_gluteus_medius_r] = 12
_name2id[_gluteus_minimus_l] = 13
_name2id[_gluteus_minimus_r] = 14
_name2id[_gracilis_l] = 15
_name2id[_gracilis_r] = 16
_name2id[_iliacus_l] = 17
_name2id[_iliacus_r] = 18
_name2id[_obturator_externus_l] = 19
_name2id[_obturator_externus_r] = 20
_name2id[_obturator_internus_l] = 21
_name2id[_obturator_internus_r] = 22
_name2id[_pectineus_l] = 23
_name2id[_pectineus_r] = 24
_name2id[_piriformis_l] = 25
_name2id[_piriformis_r] = 26
_name2id[_psoas_major_l] = 27
_name2id[_psoas_major_r] = 28
_name2id[_rectus_femoris_l] = 29
_name2id[_rectus_femoris_r] = 30
_name2id[_sartorius_l] = 31
_name2id[_sartorius_r] = 32
_name2id[_semimembranosus_l] = 33
_name2id[_semimembranosus_r] = 34
_name2id[_semitendinosus_l] = 35
_name2id[_semitendinosus_r] = 36
_name2id[_tensor_fasciae_latae_l] = 37
_name2id[_tensor_fasciae_latae_r] = 38
_name2id[_vastus_lateralis_intermedius_l] = 39
_name2id[_vastus_lateralis_intermedius_r] = 40
_name2id[_vastus_medialis_l] = 41
_name2id[_vastus_medialis_r] = 42
_name2id[_tibia_l] = 43
_name2id[_tibia_r] = 44
_name2id[_fibula_l] = 45
_name2id[_fibula_r] = 46
_name2id[_anterior_comp_l] = 47
_name2id[_anterior_comp_r] = 48
_name2id[_lateral_comp_l] = 49
_name2id[_lateral_comp_r] = 50
_name2id[_deep_posterior_comp_l] = 51
_name2id[_deep_posterior_comp_r] = 52
_name2id[_superficial_posterior_comp_l] = 53
_name2id[_superficial_posterior_comp_r] = 54
_name2id[_foot_bone_l] = 55
_name2id[_foot_bone_r] = 56
_name2id[_patella_l] = 57
_name2id[_patella_r] = 58
_name2id[_sacrum] = 59
_name2id[_skin] = 0
#
_name2id[_origin] = 60
#
_name2id[_anterior_sacral_slope] = 61
_name2id[_posterior_sacral_slope] = 62
_name2id[_rt_ASIS] = 64#63
_name2id[_lt_ASIS] = 63#64
_name2id[_rt_hip_center] = 66#65
_name2id[_lt_hip_center] = 65#66
#
_name2id[_lt_prox_femur] = 67
_name2id[_rt_prox_femur] = 68
_name2id[_lt_knee_center] = 69
_name2id[_rt_knee_center] = 70
_name2id[_lt_late_retrocondy] = 71
_name2id[_rt_late_retrocondy] = 72
_name2id[_lt_med_retrocondy] = 73
_name2id[_rt_med_retrocondy] = 74
_name2id[_lt_med_retrotibia] = 75
_name2id[_rt_med_retrotibia] = 76
_name2id[_lt_late_retrotibia] = 77
_name2id[_rt_late_retrotibia] = 78
_name2id[_lt_ankle_center] = 79
_name2id[_rt_ankle_center] = 80
_name2id[_lt_med_malleolus] = 81
_name2id[_rt_med_malleolus] = 82
_name2id[_lt_late_malleolus] = 83
_name2id[_rt_late_malleolus] = 84
#
_name2id[_talus_center_l] = 85
_name2id[_talus_center_r] = 86
#_name2id[_calcaneal_center_l] = 87
#_name2id[_calcaneal_center_r] = 88
#_name2id[_navicular_center_l] = 89
#_name2id[_navicular_center_r] = 90
#_name2id[_medial_cuneiform_center_l] = 91
#_name2id[_medial_cuneiform_center_r] = 92
#_name2id[_intermediate_cuneiform_center_l] = 93
#_name2id[_intermediate_cuneiform_center_r] = 94
#_name2id[_lateral_cuneiform_center_l] = 95
#_name2id[_lateral_cuneiform_center_r] = 96
#_name2id[_cuboid_center_l] = 97
#_name2id[_cuboid_center_r] = 98
_name2id[_talus_l] = 99
_name2id[_talus_r] = 100
#_name2id[_calcaneal_l] = 101
#_name2id[_calcaneal_r] = 102
#_name2id[_navicular_l] = 103
#_name2id[_navicular_r] = 104
#_name2id[_medial_cuneiform_l] = 105
#_name2id[_medial_cuneiform_r] = 106
#_name2id[_intermediate_cuneiform_l] = 107
#_name2id[_intermediate_cuneiform_r] = 108
#_name2id[_lateral_cuneiform_l] = 109
#_name2id[_lateral_cuneiform_r] = 110
#_name2id[_cuboid_l] = 111
#_name2id[_cuboid_r] = 112

_strcture_name_list_without_skin = [
    _sacrum,
    _pelvis_l,
    _pelvis_r,
    _femur_l,
    _femur_r,
    _adductor_l,
    _adductor_r,
    _biceps_femoris_l,
    _biceps_femoris_r,
    _gluteus_maximus_l,
    _gluteus_maximus_r,
    _gluteus_medius_l,
    _gluteus_medius_r,
    _gluteus_minimus_l,
    _gluteus_minimus_r,
    _gracilis_l,
    _gracilis_r,
    _iliacus_l,
    _iliacus_r,
    _obturator_externus_l,
    _obturator_externus_r,
    _obturator_internus_l,
    _obturator_internus_r,
    _pectineus_l,
    _pectineus_r,
    _piriformis_l,
    _piriformis_r,
    _psoas_major_l,
    _psoas_major_r,
    _rectus_femoris_l,
    _rectus_femoris_r,
    _sartorius_l,
    _sartorius_r,
    _semimembranosus_l,
    _semimembranosus_r,
    _semitendinosus_l,
    _semitendinosus_r,
    _tensor_fasciae_latae_l,
    _tensor_fasciae_latae_r,
    _vastus_lateralis_intermedius_l,
    _vastus_lateralis_intermedius_r,
    _vastus_medialis_l,
    _vastus_medialis_r,
    _tibia_l,
    _tibia_r,
    _fibula_l,
    _fibula_r,
    _anterior_comp_l,
    _anterior_comp_r,
    _lateral_comp_l,
    _lateral_comp_r,
    _deep_posterior_comp_l,
    _deep_posterior_comp_r,
    _superficial_posterior_comp_l,
    _superficial_posterior_comp_r,
    _foot_bone_l,
    _foot_bone_r,
    _patella_l,
    _patella_r
]

_strcture_name_list = copy.deepcopy(_strcture_name_list_without_skin).append(_skin)