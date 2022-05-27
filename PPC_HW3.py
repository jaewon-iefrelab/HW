import numpy as np
import pandas as pd
from scipy.stats import norm

# parameters
#T, c, h, p, gamma = 10, 1, 0.5, 10, 0.98
#mu, sigma = 20, 5

def g(y):
    s = 0
    for d in np.arange(min_,max_,step_size):
        s = s + (h*max(0,y-d)+p*max(0,d-y)) * norm.pdf(d,mu,sigma) * step_size
    return s

class BaseStockDP:
    def __init__(self):
        self.optimal_S = {t:None for t in range(T, 0, -1)}
        self.Hy = {t:{} for t in range(T, 0, -1)}
        self.y = {t:{} for t in range(T, 0, -1)}
        self.theta = {t:{} for t in range(T, 0, -1)}

    def solve(self):
        for t in range(T,0,-1):
            temp = float("inf")
            Ht = {y:None for y in range(max_, min_-1, -1)}

            if(t==T):
                for y in np.arange(min_,max_+1,1):
                    Ht[y] = (c*y + (1+gamma)*g(y))
                    if(Ht[y]<temp):
                        temp = Ht[y]
                        self.optimal_S[t] =  y

                self.Hy[t] = Ht

            else:
                for y in np.arange(min_,max_+1,1):
                    Expectation = 0
                    for d in np.arange(min_,max_,step_size):
                        theta_t = 0
                        flt_y = y-d-int(y-d)
                        int_y = int(y-d)

                        #overflow 방지
                        if(int_y>max_):
                            int_y = max_

                        if (y-d>self.optimal_S[t+1]):
                            #overflow 방지
                            if(int_y==max_):
                                theta_t = (flt_y)*self.Hy[t+1][int_y] + (1-flt_y)*self.Hy[t+1][int_y] - c*(y-d)
                            #정수가 아닌 값에 대해서는 보간법 사용하여 보정
                            else:
                                theta_t = (flt_y)*self.Hy[t+1][int_y] + (1-flt_y)*self.Hy[t+1][int_y+1] - c*(y-d)

                        else:
                            theta_t = self.Hy[t+1][self.optimal_S[t+1]] - c*(y-d)

                        Expectation = Expectation + theta_t * norm.pdf(d,mu,sigma) * step_size

                    Ht[y] = (c*y + g(y) + gamma*Expectation)
                    if(Ht[y]<temp):
                        temp = Ht[y]
                        self.optimal_S[t] =  y
                self.Hy[t] = Ht

    def answer(self):
        for t in range(T,0,-1):
            for x in list(self.Hy[t].keys()):
                self.theta[t][x] = (x<self.optimal_S[t])*(self.Hy[t][dp.optimal_S[t]]) + (x>=self.optimal_S[t])*(self.Hy[t][x]) - c*x
                self.y[t][x] = (x<self.optimal_S[t])*(self.optimal_S[t]) + (x>=self.optimal_S[t])*x

def main():
    input_list = input("T, c, h, p, gamma ,mu, sigma = ")
    param_list = [float(i) for i in input_list.split(',')]
    T = param_list[0]
    C = param_list[1]
    h = param_list[2]
    p = param_list[3]
    gamma = param_list[4]
    mu = param_list[5]
    sigma = param_list[6]
    max_ = mu+6*sigma
    min_ = mu-6*sigma
    step_size = 0.001*sigma
    print("Solve Start")
    dp = BaseStockDP()
    dp.solve()
    dp.answer()
    print("Solve End")

    save = input("Want to save?(Y/N)")
    if save == 'Y':
            pd.DataFrame([dp.optimal_S]).to_csv('optimal_S.csv',index=True)
            pd.DataFrame(dp.theta).to_csv('theta.csv',index=True)
            pd.DataFrame(dp.y).to_csv('y.csv',index=True)


if __name__ == '__main__':
    main()
