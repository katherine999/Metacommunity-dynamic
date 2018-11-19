# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 16:41:01 2018

@author: Lin
"""
import numpy as np
import random
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm

def generate_habitat(e, length, width):
    microsite_e_values = np.random.normal(loc=0, scale=0.025, size=(length, width))+e
    microsite_individuals = [[None for i in range(length)] for i in range(width)]
    return {'microsite_e_values':microsite_e_values, 'microsite_individuals':microsite_individuals}
# habitat里包括了microsite_e_values 和 microsite_individuals 两个值, 分别是5*5的list
   
def generate_patch():
    h0 = generate_habitat(e=0.2, length=5, width=5)
    h1 = generate_habitat(e=0.4, length=5, width=5)
    h2 = generate_habitat(e=0.6, length=5, width=5)
    h3 = generate_habitat(e=0.8, length=5, width=5)
    patch = {'h0':h0, 'h1':h1, 'h2':h2, 'h3':h3}
    return patch
# patch 里包含了四个habitat的集合，h0、h1、h2、h3
   
def generate_matacommunity():
    metacommunity = {}
    for i in range(9):
        patch = generate_patch()
        metacommunity['patch%s'%(str(i))] = patch
    return metacommunity
# matacommunity里包含了9个patch，p1——p9
    
def init_individual(e, L=10, sexual=None):               # sexual为None, male, female
    phenotype = e + random.gauss(0,0.025)
    genotype = [0 if e<np.random.uniform(0,1,1) else 1 
                          for i in range(L)]
    return {'identifier': int(e/0.2), 'sexual': sexual, 'phenotype':phenotype, 'genotype':genotype}
# 根据e值，初始化一个个体的phenotype和genotype值，表示该初始化物种适应对应的habitat环境
# 一个个体包含以下信息， identifier表示 不同的物种表识，sexual表示性别（有性繁殖分为male和female，无性繁殖则为None）
# phenotype为表型，genotype为基因型    
   
    
def mutation(μ, individual):                             # 突变的函数
    genotype = individual['genotype']
    phenotype = individual['phenotype']
    for allelic in range(len(genotype)):
        if μ > np.random.uniform(0,1,1):                 # 每个基因以μ的概率突变成它的一个等位基因
            if genotype[allelic] ==0: 
                genotype[allelic] = 1
                #print('mutation to 1')
            if genotype[allelic] ==1: 
                genotype[allelic] = 0 
                #print('mutation to 0')
    phenotype = np.mean(genotype) + random.gauss(0,0.025)   # 表现型由基因和一些非遗传因素决定
    individual['genotype']=genotype
    individual['phenotype']=phenotype
    return individual

def parents():
    pass 
def reproduce(μ, parents, reproduction_type):    #当无性繁殖时parents为一个亲本，当有性繁殖是parents为两个亲本
    if reproduction_type == 'asexual':           # 无性繁殖
        new = mutation(μ, parents)               # 突变发生在繁殖后代的过程中
         
    if reproduction_type == 'sexual':            # 有性繁殖
        female = parents[0]
        male = parents[1]
        m = random.sample([i for i in range(10)],5)  
        n_identifier = female['identifier']
        n_genotype = []
            
        if 0.5 > np.random.uniform(0,1,1):       # 后代的性别，female与male均为50%
            n_sexual = 'female'
        else:
            n_sexual = 'male'

        for allelic in range(len(female['genotype'])): # 5个基因来自母本，5个基因来自父本
            if allelic in m:
                n_genotype.append(male['genotype'][allelic])
            else:
                n_genotype.append(female['genotype'][allelic])
                
        n_phenotype = np.mean(n_genotype) + random.gauss(0,0.025)
        new = {'identifier':n_identifier, 'sexual':n_sexual,'phenotype':n_phenotype, 'genotype':n_genotype}
        new = mutation(μ, new)                    # 突变发生在繁殖后代的过程中
    return new

def survival_rate(d, z, e, w = 0.5):
    survival_rate = (1-d) * math.exp((-1)*math.pow(((z-e)/w),2))       # 存活的几率符合正态分布
    print(survival_rate)
    return survival_rate
# d为基本的死亡量，z为某个个体的表型，e为该microsite的环境值，w为width of the fitness function
    
def find_pools(patch_id, immigrant_pool):  # patch_id 当前的patch，即dispersal within的patch
    among = []
    for p_id in immigrant_pool:
        I_p = immigrant_pool[p_id]
        if p_id == patch_id:               # 当前遍历的patch为within_patch
            within = I_p              # 此patch的后代分类为within_pool
        else:
            among += I_p              # 其他patch的后代分类为among_pool
    return within, among  # 返回dispersal within 和 amomg的新个体

######################################     模拟算法     ###########################################




#####################################    初始化    #######################################################################
def initialization():
    for patch_id in metacommunity:       # metacommunity为一个字典，此语法得到的patch_id，为字典的key值
        if int(patch_id[5:]) >= init_occupied_patch_num: # 模拟开始时，已被占有的patches的数量
            break
        else:
            patch = metacommunity[patch_id] 
            for habitat_id in patch:     # patch为一个字典，此语法得到habitat_id, 为字典的一个key值
                '''habitat = patch[habitat_id]
                microsite_individuals = habitat['microsite_individuals'] # 一个5*5的list'''
                e = (int(habitat_id[1:]) + 1) * 0.2
                for length in range(5):
                    for width in range(5):
                        new_individual = init_individual(e, L=10, sexual=None)
                        # e值为该habitat的e值，L为10对等位基因，繁殖方式为无性繁殖，individual没有性别None
                        metacommunity[patch_id][habitat_id]['microsite_individuals'][length][width] = new_individual
 
                    
#####################################     死亡选择     #####################################################################
def eco_evolutionary_dynamic():
    for patch_id in metacommunity:      # metacommunity为一个字典，此语法得到的patch_id，为字典的key值
        patch = metacommunity[patch_id]
        for habitat_id in patch:        # patch为一个字典，此语法得到habitat_id, 为字典的一个key值
            habitat = patch[habitat_id]
            microsite_individuals_set = habitat['microsite_individuals']   # 一个5*5的list，包含25个e值
            microsite_e_values_set = habitat['microsite_e_values']         # 一个5*5的list，包含25个individual值
            for length in range(5):
                for width in range(5):
                    microsite_e_value = microsite_e_values_set[length][width]  # 一个microsite的e_value
                    if microsite_individuals_set[length][width] != None:
                        microsite_individual = microsite_individuals_set[length][width] #一个microsite的individual值
                        phenotype = microsite_individual['phenotype']                  # 该microsite个体的一个phenotype值
                        survival_p = survival_rate(d=0.1, z=phenotype, e=microsite_e_value, w = 0.5)
                        # 通过表型和环境的e_value计算该个体的存活率
                        if survival_p > np.random.uniform(0,1,1):
                            metacommunity[patch_id][habitat_id]['microsite_individuals'][length][width] = None
                            # 表示该个体已经死亡，用None表示
                
                
    ##################################     繁殖后代     ######################################################################                

    offsprings_pool = {}      

    # 无性繁殖
    birth_rate = 0.5                          # 无性繁殖时，出生率为0.5

    for patch_id in metacommunity:            # metacommunity为一个字典，此语法得到的patch_id，为字典的key值
        patch = metacommunity[patch_id]
        pool = []                             # 每个patch里分别有一个offspring的pool
        for habitat_id in patch:              # patch为一个字典，此语法得到habitat_id, 为字典的一个key值
            habitat = patch[habitat_id]

            microsite_individuals_set = habitat['microsite_individuals'] # 一个5*5的list，包含25个individual值
            #print('habitat_id', habitat_id)
            for length in range(5):
                for width in range(5):
                    if birth_rate < np.random.uniform(0,1,1) and microsite_individuals_set[length][width] != None:
                        # 该microsite不为empty状态，并且以某一birth rate繁殖后代
                        parent = microsite_individuals_set[length][width]
                        new_individual = reproduce(μ=0.0001, parents=parent, reproduction_type='asexual')
                        #print('new_individual=', new_individual['genotype'], new_individual['phenotype'])
                        pool.append(new_individual)          # 将一个patch里的所有后代储存起来
                        offsprings_pool[patch_id] = pool     # 每个patch的子代个体集合
                                                         # 并以patch_id作为该集合的一个key值
        if pool==[]:
            offsprings_pool[patch_id] = [None]                # 用None表示该patch中个体无子代个体
                    
                    
    ####################################     个体迁移      ####################################################################333

    immigrant_pool = offsprings_pool
    m_within = 0.1
    m_among = 0.01



    for patch_id in metacommunity:                  # metacommunity为一个字典，此语法得到的patch_id，为字典的key值                                                          
        patch = metacommunity[patch_id]
        within_pool, among_pool = find_pools(patch_id, immigrant_pool)     
        # 将所有新的后代分为within、among两大类
        for habitat_id in patch:                    # patch为一个字典，此语法得到habitat_id, 为字典的一个key值
            habitat = patch[habitat_id]
            microsite_individuals_set = habitat['microsite_individuals']  # 一个5*5的list，包含25个individual值
            for length in range(5):
                for width in range(5):
                    if microsite_individuals_set[length][width] == None: # 找出当前处于empty状态的patch    
                        if m_within < np.random.uniform(0,1,1) and within_pool != None:          # m_within的概率在同一个patch之间迁移
                            dispersal = random.sample(within_pool, 1)[0]
            
                            metacommunity[patch_id][habitat_id]['microsite_individuals'][length][width] = dispersal
                        
                        elif m_among < np.random.uniform(0,1,1) and among_pool != None:         # m_among的概率在不同的patch之间迁移
                            dispersal = random.sample(among_pool, 1)[0]

                            metacommunity[patch_id][habitat_id]['microsite_individuals'][length][width] = dispersal
####################################### 可视化 ##################################################################
def sub_hotpower(z_set):
    x, y = np.random.rand(25), np.random.rand(4)
    plt.imshow(z_set, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),cmap=cm.hot, norm=LogNorm())    
    
def get_z_set(patch):
    z = []
    for habitat_id in patch:
        habitat = patch[habitat_id]
        for length in range(5):
            for width in range(5):
                if habitat['microsite_individuals'][length][width]!=None:
                    z.append(habitat['microsite_individuals'][length][width]['phenotype'])
                else:
                    z.append(0.001)
    z = np.array(z).reshape(4,25)
    return np.abs(z)

def show_phenotype():
    fig = plt.figure()
    fig.add_subplot(331)
    patch0 = metacommunity['patch0']
    z_set = get_z_set(patch0)
    sub_hotpower(z_set)

    fig.add_subplot(332)
    patch1 = metacommunity['patch1']
    z_set = get_z_set(patch1)
    sub_hotpower(z_set)

    fig.add_subplot(333)
    patch2 = metacommunity['patch2']
    z_set = get_z_set(patch2)
    sub_hotpower(z_set)

    fig.add_subplot(334) 
    patch3 = metacommunity['patch3']
    z_set = get_z_set(patch3)
    sub_hotpower(z_set)

    fig.add_subplot(335)
    patch4 = metacommunity['patch4']
    z_set = get_z_set(patch4)
    sub_hotpower(z_set)

    fig.add_subplot(336)
    patch5 = metacommunity['patch5']
    z_set = get_z_set(patch5)
    sub_hotpower(z_set)

    fig.add_subplot(337)
    patch6 = metacommunity['patch6']
    z_set = get_z_set(patch6)
    sub_hotpower(z_set)

    fig.add_subplot(338)
    patch7 = metacommunity['patch7']
    z_set = get_z_set(patch7)
    sub_hotpower(z_set)

    fig.add_subplot(339) 
    patch8 = metacommunity['patch8']
    z_set = get_z_set(patch8)
    sub_hotpower(z_set)
    
    plt.colorbar()
    plt.show()    


#######################################     运行模型     #########################################################################

metacommunity = generate_matacommunity()
# metacommunity里包含了9个patches, 每个patch里包含立4个habitat
# habiat里储存了两个信息，microsit_e_value 和 individuals的值
init_occupied_patch_num = 8
# 一开始occupied的patch的数量，可以是1或8
initialization()
show_phenotype()


for time_step in range(1):
    eco_evolutionary_dynamic() 
    if time_step % 2 == 0:
        print('time_step=', time_step+1)
        show_phenotype()
                    

    
   
    
    
    




    
    
    
    
    