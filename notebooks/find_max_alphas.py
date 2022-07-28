
import numpy as np
from classy import Class

#Try computing at alpha1, alpha2. 
#If both don't work, set new alpha1=alpha1-0.1, and new alpha2= old alpha1. 
#If alpha2 doesn't work but alpha1 does: Set alpha2 = (alpha2+alpha1)/2 and try again. 
#If alpha1 and alpha2 work, make new alpha1 = alpha2, and new alpha2 = 2 * alpha1 - alpha2. 


def test_alpha(r,deltaN,mp,me,alpha):
    try:
        settings = {'output':'mPk, tCl',
           'omega_b':0.0224,
           'omega_cdm':0.119,
           '100*theta_s':1.04,
           'ln10^{10}A_s':3.05,
           'n_s':0.965,
           'tau_reio':0.0576,
            'YHe_twin':0,
           'r_all_twin':r,
           'Delta_N_twin':deltaN,
           'm_p_dark':mp,
            'm_e_dark':me,
            'alpha_dark':alpha}
        M = Class()
        M.set(settings)
        M.compute()
    except:
        return False
    else:
        return True

def find_max_alpha(r,deltaN,mp,me):

    alpha_guess_1 = 0.05
    alpha_guess_2 = 0.06

    
    good_1 = test_alpha(r,deltaN,mp,me,alpha_guess_1)
    
    good_2 = test_alpha(r,deltaN,mp,me,alpha_guess_2)
    

    while not (((alpha_guess_2 - alpha_guess_1)/(alpha_guess_2) < 0.001 ) and (good_1 and (not good_2))):
        if (not good_1) and (not good_2):
            alpha_guess_2 = alpha_guess_1
            alpha_guess_1 = alpha_guess_1 / 2

        elif (good_1) and (not good_2):
            alpha_guess_2 = (alpha_guess_2 + alpha_guess_1) / 2

        elif (not good_1) and (good_2):
            #print('This should not be able to happen, guess_1 should be < guess_2')
            #print(alpha_guess_1)
            #print(alpha_guess_2)
            return alpha_guess_1+100
        elif good_1 and good_2: 
            alpha_guess_1_new = alpha_guess_2
            alpha_guess_2_new = 2 * alpha_guess_2 - alpha_guess_1
            alpha_guess_1 = alpha_guess_1_new
            alpha_guess_2 = alpha_guess_2_new

        good_1 = test_alpha(r,deltaN,mp,me,alpha_guess_1)

        good_2 = test_alpha(r,deltaN,mp,me,alpha_guess_2)

    return alpha_guess_1

#params_list = [[0.1,0.1,0.1,0.0001],[0.1,0.1,0.1,0.001],[0.1,0.1,0.1,0.01],[0.1,0.1,10,0.0001],[0.1,0.1,10,0.001],[0.1,0.1,10,0.01],[0.1,0.01,1,0.001],[0.1,0.3,1,0.001],[0.1,0.1,1,0.0001],[0.1,0.1,1,0.001],[0.1,0.1,1,0.01],[0.1,0.1,1,0.05]]
params_list = [[0.1,0.1,1,0.0001],[0.1,0.1,1,0.0005],[0.1,0.1,1,0.0007],[0.1,0.1,1,0.0008],[0.1,0.1,1,0.001],[0.1,0.1,1,0.002],[0.1,0.1,1,.003],[0.1,0.1,1,0.004]]
f=open('/project/d/dcurtin/jpbarron/ADM/highest_alphas.txt','w')
for params in params_list:
    alpha=find_max_alpha(*params)
    params.append(alpha)
    f.write(','.join([str(x) for x in params]))
    
f.close()
    

