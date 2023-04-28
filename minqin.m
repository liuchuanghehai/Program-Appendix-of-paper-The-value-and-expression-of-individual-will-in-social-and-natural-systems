clear
clc
%% 数据输入
I = [];
I(1) = 102.8858;
Q = [5.1
4.6
6.2
5.7
4.49
4.88
3.7
3.74
5.1
3.5
4
6
4.9
4.26
4.45
4.28
3.47
3.4
3.16
2.95
2.84
2.77
2.4
2.18
2.38
2.45
2.18
2.6
2.24
2.15
2.29
2.22
2.1
2.28
1.77
1.57
1.21
2.26
1.72
1.51
1.55
1.41
1.07
1.04
1.15
1.18
0.51];
P=[65.36
106.06
87.53
87.53
28.72
12.14
63.09
76.49
39.14
76.34
64.92
126.91
77.82
109.65
88.89
59.47
163.27
39.12
72.2
34.63
106.13
140.43
79.93
68.45
35.18
23.84
87.5
49.8
114.44
39.03
82.51
79.69
23.82
26.52
49.26
59.77
87.77
158.36
68.89
94.74
56.17
65.3
59.46
84.22
55.96
119.07
78.75];
Pop=[0.02
-0.06633811
-0.072262489
-0.072262489
-0.072262489
0.016279755
0.016279755
0.016279755
0.016279755
0.030593352
0.030593352
0.030593352
0.030593352
0.030593352
0.018645514
0.018645514
0.018645514
0.018645514
0.018645514
0.002745412
0.002745412
0.002745412
0.002745412
0.002745412
0.006250416
0.006250416
0.004486768
0.004486768
0.004486768
0.00590065
0.016553915
0.012742718
0.013048399
0.017294899
0.010852589
0.015111703
0.00697677
0.015184446
0.010374857
0.019539933
0.032978626
0.015510104
0.009386524
0.006267021
0.006194985
0.006803033
0.000661371];
years=[];
years(1) = 1956;
ki = [];
S = [];
S(1) = Q(1)*10000 + 2/3*P(1)*150 - 1.1*(400*60/62+200*2/60)*I(1);
wpp = [];
wpp(1) = 2;
po = 24;
c = 200;
record=[];
wmax=[];
dss = [];
%% 计算
%for i = 2:64 
    i = 2;%1957年到2020年
while i+1957<=2050
    years(i) = i+1956;%年份迭代
   %% 灌溉水利用系数
    if i+1957<1965
        n = 0.33+0.01*i;
    elseif and(i+1957>1965,i+1957<2005)   
        n = 0.54;
    elseif and(i+1957>=2005,i+1957<2020)
        n = n+0.012;
        if n>0.6
            n = 0.6;
        end
    end
    %% 社会扩张意愿
    if i+1957<1980
        faI = 0.6387;
    elseif i+1957>=1980
        faI = 0.6387-0.03746*(i+1957-1980);
        if faI<=0.4514
            faI = 0.4514;
        end
    end
    %% 种植结构
    if i+1957<1980
        Ig = 60/62;
    elseif and(i+1957>=1980,i+1957<1992)
        Ig = 60/62-0.012*(i+1957-1980+1);
    else
        Ig = 43/53-0.0234*(i+1957-1992+1);
        if Ig<=0.4
            Ig = 0.4;
        end
    end
    %% 最大增长率
    if and(i+1957>1978,i+1957<=1990)
        r = 0.05;
    elseif i+1957>1990
        r = 0.08;
    elseif i+1957<=1978
        r = 0.08;
        if and(i+1957>1961,i+1957<1966)
        r = 0.02;
        end
    end
    if i+1957<2005
        mg = 480;
        me = 250;
    else
        mg = mg-46;
        if mg <=250
            mg = 250;
        end
    end
    %% 地下水
    if and(i+1957>1965,i+1957<=2005)
        if wp <6
            wp = wp + min(0.4,wp*(1-wp/6)*((I(i-1))*(mg*Ig+me*(1-Ig))/(n*Q(i-1)*10000+0.6*10000*wp+0.6*2/30000*P(i-1)*(I(i-1)))));
            %wp = wp + min(0.4,max((0.02*150*faI+I(i-1))*(550*Ig+250*(1-Ig))-(n*Q(i)*10000+0.8*2/3*P(i)*(0.02*150*faI+I(i-1))),0));
            if wp >6
                wp =6;
            end
        else
            wp = 6;
        end
    elseif i+1957<=1965 
        wp = 2;
    elseif and(i+1957>2005,i+1957<2010)
        wp = wp - 0.514;
    elseif i+1957>=2010
        wp = 0.89;
    end
    wpp(i) = wp;
    %% 降雨与来水
    if and(i+1957>2003,i+1957<2006)
        Q(i) = Q(i-1);
        P(i) = 73+randn(1)*30;
    elseif and(i+1957>2005,i+1957<2010)
        Q(i) = 0.5+0.24*(i+1957-2005);
        P(i) = 73+randn(1)*30;
        if Q(i)>=2.9
            Q(i) = 2.9+randn(1)*0.5 ;
        end
    elseif i+1957>=2010
        Q(i) = 2.9+randn(1)*0.5 ;
        P(i) = 73+randn(1)*30;
    end
    %% 耕地扩张速度测算
    dI = I(i-1)*min(r*(1.2-I(i-1)/c/faI),r*(1+Pop(i-1)-I(i-1)/c)*(1-1.1*(mg*Ig+me*(1-Ig))*I(i-1)/((n*Q(i)*10000+0.6*wp*10000+0.6*2/30000*P(i)*I(i-1)))));
    c = c*(1-0.003);%%土地盐碱率
    if i+1957>2010
        c = c/(1-0.003);
    end
    %1.2-I(i-1)/150/faI   *(2*po/Ig/I(i-1))
    %% 人口
    if i+1957<2003
        po = po*(Pop(i)+1);%%人口
    elseif i+1957>=2010
        if po>25;
            Pop(i) = 0.003-0.01;
            ui = 0.01;
        else
            ui = 0;
            Pop(i) = 0.003;
        end
        po = po*(1+0.003-ui);
    else 
        Pop(i) = 0.003;
        ui = 0;
        po = po*(1+0.003-ui);
    end
    ki(i) = dI/I(i-1);
    %*(0.001*Pop(i-1)-0.01*UI)  *(0.001*Pop(i)) 0.0013*
    %*1.4*po/(Ig*I(i-1))
    %% 水量
    I(i) = I(i-1)+dI;
    ds = -(n*Q(i)*10000 + 0.6*2/30000*P(i)*I(i))/0.6 + 1.1*(mg*Ig+me*(1-Ig))*I(i)/0.6;
    wmax(i) = min(ds,wp*10000);
    S(i) = ds;
    dss(i) = Q(i)*10000+2/30000*P(i)*I(i)-min(n*Q(i)*10000+0.6*2/30000*P(i)*I(i)+wmax(i),1.1*(mg*Ig+me*(1-Ig))*I(i));
    i = i+1;
end
plot(years,I)