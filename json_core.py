import json
def config(arr,arr2,arr3,arr4,arr5,arbarr,sarbarr,darr,trarr,oparr,filename):
    data = {
        "driver_settings" : {
            "Ibias" : float(arr[0])*1e-12,
            "Iptp" : float(arr[1])*1e-12,
            "number of Pulses" : arr[2],
            "integration rate" : arr[3],
            "filter bandwidth" : arr[4],
            "filter order" : arr[5],
            "pulse repetiotion rate" : arr[6],
            "pulse width" : arr[7],
            "decay time" : arr[8]
        }   , 
        "laser_settings" : {
            "eta" : arr[9] , #0.3,
            "t_ph" : arr[10],
            "N_transparency" : arr[11],
            "N_threshhold" : arr[12],
            "chi" : arr[13],
            "C_sp" : arr[14],
            "Gamma" : arr[15],
            "tau_e" : arr[16],
            "alpha" : arr[17],
            "rh" : arr[18],
            "th" : arr[19],
            "lambda" : arr[20],
            "kap_om" : arr[21],
            "t1" : arr[22],
            "w_M" : arr[23],
            "om_th" : arr[24],
            "Ll" : arr[25],
            "kap_inj" : arr[26],
            "r1" : arr[27],
            "r2" : arr[28],
            "rhs" : arr[29],
            "k12" : arr[30],
            "tauj" : arr[31],
            "lambdas": arr[32],
            "Cj" : arr[33],
            "int_delay" : arr[34],
            "phi1" : arr[35] ,
            "phi2" : arr[36]            
        }   ,
        "wdm_settings" : {
            "filter bandwidth" : arr2[0],
            "filter order" : arr2[1],
            "filter shift" : arr2[2]
        } ,
        "interferometer settings": {
            "int_delay" : arr[34],
            "phi1" : arr3[35] ,
            "phi2" : arr3[36]   
        },
        "phot_settings" : {
            "phot bandwidth" : arr5[0],
            "phot order" : arr5[1],
            "phot shift" : arr5[2]
        },
        "slave_driver_settings" : {
            "Ibias" : float(arr3[0])*1e-12,
            "Iptp" : float(arr3[1])*1e-12,
            "number of Pulses" : arr3[2],
            "integration rate" : arr3[3],
            "filter bandwidth" : arr3[4],
            "filter order" : arr3[5],
            "pulse repetiotion rate" : arr3[6],
            "pulse width" : arr3[7],
            "decay time" : arr3[8]
        }   ,
        "slave_laser_settings" : {
            "eta" : arr3[9] , #0.3,
            "t_ph" : arr3[10],
            "N_transparency" : arr3[11],
            "N_threshhold" : arr3[12],
            "chi" : arr3[13],
            "C_sp" : arr3[14],
            "Gamma" : arr3[15],
            "tau_e" : arr3[16],
            "alpha" : arr3[17],
            "rh" : arr3[18],
            "th" : arr3[19],
            "lambda" : arr3[20],
            "kap_om" : arr3[21],
            "t1" : arr3[22],
            "w_M" : arr3[23],
            "om_th" : arr3[24],
            "Ll" : arr3[25],
            "kap_inj" : arr3[26],
            "r1" : arr3[27],
            "r2" : arr3[28],
            "rhs" : arr3[29],
            "k12" : arr3[30],
            "tauj" : arr3[31],
            "lambdas": arr3[32],
            "Cj" : arr3[33],
            "int_delay" : arr3[34],
            "phi1" : arr3[35] ,
            "phi2" : arr3[36]          
        }, 
        "slave_wdm_settings" : {
            "filter bandwidth" : arr4[0],
            "filter order" : arr4[1],
            "filter shift" : arr4[2]
        },
        "arbpulse_settings":{ 
            "Ibias" : arbarr[0], 
            "Iptp": arbarr[1], 
            "pulse_width": arbarr[2], 
            "gap_current" : arbarr[3], 
            "gap_width" : arbarr[4], 
            "gap_position" : arbarr[5], 
            "gap_repetition_rate" : arbarr[6], 
            "time_shift" : arbarr[7], 
            "integration_rate": arbarr[8], 
            "incline" : arbarr[9]
        },
        "slavearbpulse_settings":{
            "Ibias" : sarbarr[0], 
            "Iptp": sarbarr[1], 
            "pulse_width": sarbarr[2], 
            "gap_current" : sarbarr[3], 
            "gap_width" : sarbarr[4], 
            "gap_position" : sarbarr[5], 
            "gap_repetition_rate" : sarbarr[6], 
            "time_shift" : sarbarr[7], 
            "integration_rate": sarbarr[8], 
            "incline" : sarbarr[9]
        },
        "diagram_settings":{
            "Qinj_min" : darr[0],
            "Qinj_max" : darr[1],
            "phi_min" : darr[2],
            "phi_max" : darr[3],
            "sol_point" : darr[4],
            "pointsnum" : darr[5]
        },
        "transfer_settings":{
             "freq_range_start" : trarr[0],
             "freq_range_end" : trarr[1]
        },
        "option_settings":{
            "driver_filter" : oparr[0],
            "wdm_filter" : oparr[1],
            "wdm_s_filter" : oparr[2],
            "calc_scheme" : oparr[3]
        }
    }

    return data
    #with open(filename,"w") as write_file:
    #    json.dump(data,write_file)
        #json.dump(laser_settings,write_file)
def conf_pulsetrain(pulseid,slavepulseid,array,slavearray,filename = None):
    data = {
        "pulse_num" : pulseid,
        "slavepulse_num" : slavepulseid
    }
    train = []
    slavetrain = []
    for i in range(0,pulseid):
        pulsedata = {
            "arbpulse_settings":{ 
                "Ibias" : array[0][i], 
                "Iptp": array[1][i], 
                "pulse_width": array[2][i], 
                "gap_current" : array[3][i], 
                "gap_width" : array[4][i], 
                "gap_position" : array[5][i], 
                "gap_repetition_rate" : array[6][i], 
                "time_shift" : array[7][i], 
                "integration_rate": array[8][i], 
                "incline" : array[9][i]
            }
        }
        train.append(pulsedata)
    for i in range(0,slavepulseid):
        pulsedata = {
            "slavearbpulse_settings":{ 
                "Ibias" : slavearray[0][i], 
                "Iptp": slavearray[1][i], 
                "pulse_width": slavearray[2][i], 
                "gap_current" : slavearray[3][i], 
                "gap_width" : slavearray[4][i], 
                "gap_position" : slavearray[5][i], 
                "gap_repetition_rate" : slavearray[6][i], 
                "time_shift" : slavearray[7][i], 
                "integration_rate": slavearray[8][i], 
                "incline" : slavearray[9][i]
            }
        }
        slavetrain.append(pulsedata)
    return [data,train,slavetrain]