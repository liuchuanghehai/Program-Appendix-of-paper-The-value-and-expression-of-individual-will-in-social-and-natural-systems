clear
clc
%% 数据输入
I = [];
I(1) = 131.8005;
Q = [9.435 
8.465 
12.604 
10.302 
8.227 
10.216 
6.159 
6.845 
10.688 
5.829 
6.455 
11.367 
8.029 
8.057 
9.343 
7.566 
6.893 
7.887 
6.820 
7.753 
8.956 
8.862 
7.563 
8.186 
7.994 
9.158 
7.526 
9.536 
7.748 
7.518 
7.666 
6.700 
10.291 
10.843 
8.231 
5.166 
7.209 
10.230 
7.222 
7.096 
6.948 
7.064 
6.884 
6.002 
7.870 
6.061 
6.499 
11.243 
7.790 
8.492 
9.360 
9.894];
P=[19.03
77
75.55
29.24
125.1
96.9
24.9
49.3
132.88
94.54
131.51
139.8
161.5
96.75
163
119.2
52.3
114.1
49.8
74.5
52.8
107.2
66.8
100.68
79.2
43.2
68
136.53
87.8
180.4
88.9
79.6
110.7
71.1
102.7
38.5
100.7
175.6
120.3
151.9
147.8
117.5
87.6
121.2
129
116.9
117
99.8
100
90
110
120
100];
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
S(1) = Q(1)*10000 + 2/3*P(1)*150 - 1*(400*60/62+200*2/60)*I(1);
wpp = [];
wpp(1) = 1;
po = 24;
c = 300;
record=zeros(9,4);
wmax=[];
dss = [];
qout = [];
qout(1) = 5.1;
%% 计算
%for i = 2:64 
    i = 2;%1957年到2020年
    n = 0.5;
    faqout = 1;
    ui = 0;
    elsewater =1;
while i+1957<=2020
    years(i) = i+1956;%年份迭代
   %% 敏感性部分
    if i+1956==2008
        n_s = n;
        r1_s = r;
        r2_s = 0.4;
        %m_s = mg*Ig+me*(1-Ig);
        ig_s = Ig;
        mg_s = mg;
        me_s = me;
        c_s = c;
        wmax = 4.8;
        if i+1957>2005
            wmax = 4.18;
        end
        para = [n_s,r1_s,r2_s,ig_s,mg_s,me_s,c_s,wmax];
        y = [0.5,0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.5];
        check = zeros(7,8);
        for index = 1:10
        %index=fix(10*rand(1))+1;
        check(1,:) = para;
        result = zeros(9,4);
        for j = 2:9
            check(j,:) = check(j-1,:);
            check(j,j-1) = para(j-1)*y(index);
        end
        for j = 1:9
            n_s = check(j,1);
            r1_s = check(j,2);
            r2_s = check(j,3);
            %m_s = check(j,4);
            ig_s = check(j,4);
            mg_s = check(j,5);
            me_s = check(j,6);
            c_s = check(j,7);
            wmax = check(j,8);
            %% 地下水最大抽取量
                if wp <wmax
                    if Q(i)*10000+0.8*P(i-1)*I(i-1)*2/3/n_s+0.5*wp/n_s-1.1*(mg_s*ig_s+(1-ig_s)*me_s)*(I(i-1)*(1+r1_s*(1.2-I(i-1)/c_s/faI)))/n_s>=0
                        wp = wp+0;
                    else
                        wp = wp + min(r2_s,wp*(1-wp/wmax+((I(i-1)*(1+r*(1.2-I(i-1)/c_s/faI)))*1*(mg_s*ig_s+(1-ig_s)*me_s)/(n_s*Q(i-1)*10000+0.6*10000*wp+0.8*P(i-1)*I(i-1)*2/3))));
                    end
                    if wp >=wmax
                        wp =wmax;
                    end
                else
                    wp = wmax;
                end
            if i+1957<=1960 
                wp = 1;
            elseif and(i+1957>2005,i+1957<2010)
                wp = wp - (wp-4.18)/5;
            elseif i+1957>=2010
                wp = 4.18;
            end
            wpp(i) = wp;
            %% 灌溉面积模拟
            dI = I(i-1)*min(r1_s*(1.2-I(i-1)/c_s/faI),r*(1-I(i-1)/c_s)*(1.2-1.1*(mg_s*ig_s+(1-ig_s)*me_s)*I(i-1)/((n_s*(Q(i)-0.21*faqout*Q(i))*10000+0.5*wp*10000+0.8*P(i-1)*I(i-1)*2/3))));
            c = c*(1-0.003);%%土地盐碱率+Pop(i-1)
            I(i) = I(i-1)+dI;
            %% 实际用水量
            costq(i) = (1.1*mg_s*ig_s+(1-ig_s)*me_s)*I(i);
            ds = (1.1*mg_s*ig_s+(1-ig_s)*me_s)*I(i)/n_s-1*0.5*wp*10000/n_s-0.8*P(i-1)*I(i-1)*2/3/n_s;%这里0.8是纯井观区的取水，非纯灌区优先地下水
            S(1) = 8.825*10000;
            downwater(i) = wp*(1);
            S(i) = elsewater*ds;
            dss(i) = max((Q(i)+0.21*faqout*Q(i))*10000-S(i),0);
            result(j,1) = dI;
            result(j,2) = costq(i)-costq(i-1);
            result(j,3) = S(i)-S(i-1);
            result(j,4) = downwater(i)-downwater(i-1);
        end
        record = abs(result/abs(y(index)-1))+record;
        end
        o = i*9-12;
        xlswrite('name5.xlsx',record/11,'Sheet1',num2str(o));         
    %% 非敏感性部分
    else
   %% 灌溉水利用系数
    if i+1957<2005
        n = 0.5;
    elseif i+1957>2005
        n = n+0.0067;
        if n>0.6
            n = 0.6;
        end
    end
    %% 社会扩张意愿
    if i+1957<1966
        faI = 0.6387;
    elseif and(i+1957>=1966,i+1957<=1976)
        faI = 0.38;
    elseif i+1957>1976
        faI = faI+0.01;
        if faI>=0.6387
            faI = 0.6387;
        end
    end
    %% 种植结构
    if i+1957<1980
        Ig = 1;
    elseif and(i+1957>=1980,i+1957<2005)
        Ig = Ig-0.01;
    else
        Ig = Ig-0.017;
        if Ig<=0.5
            Ig = 0.5;
        end
    end
    %% 最大增长率
    if and(i+1957>1978,i+1957<=1990)
        r = 0.05;
    elseif i+1957>1990
        r = 0.08;
    elseif i+1957<=1978
        r = 0.08;
%         if and(i+1957>1961,i+1957<1966)
%         r = 0.02;
%         end
    end
    if i+1957<2005
        mg = 300;
        me = 250;
    end
    %% 地下水
    if and(i+1957>=1960,i+1957<=2000)
        if wp <5
            if Q(i)*10000+0.8*P(i-1)*I(i-1)*2/3/n+0.5*wp/n-1.1*((mg*Ig+me*(1-Ig)))*(I(i-1)*(1+r*(1.2-I(i-1)/c/faI)))/n>=0
                wp = wp+0;
            else
                wp = wp + min(0.4,wp*(1-wp/4.8+((I(i-1)*(1+r*(1.2-I(i-1)/c/faI)))*1*((mg*Ig+me*(1-Ig)))/(n*Q(i-1)*10000+0.6*10000*wp+0.8*P(i-1)*I(i-1)*2/3))));
                %wp = wp + min(0.4,max((0.02*150*faI+I(i-1))*(550*Ig+250*(1-Ig))-(n*Q(i)*10000+0.8*2/3*P(i)*(0.02*150*faI+I(i-1))),0));
            end
            if wp >=4.8
                wp =4.8;
            end
        else
            wp = 4.8;
        end
    elseif i+1957<=1960 
        wp = 1;
    elseif and(i+1957>2005,i+1957<2010)
        wp = wp - (wp-4.18)/5;
    elseif i+1957>=2010
        wp = 4.18;
    end
    wpp(i) = wp;
    %% 降雨与来水
    if i+1957>2007 
        Q(i) = 8+randn(1)*1.58 ;
        P(i) = 98+randn(1)*39.58;
    end
    %% 耕地扩张速度测算
    if i+1957<1980
        faqout = faqout-0.1;
        if faqout <=0
            faqout = 0;
        end
    elseif i+1957>2005
        faqout = faqout + 0.3;
        if faqout>=1
            faqout = 1;
        end
    end
    dI = I(i-1)*min(r*(1.2-I(i-1)/c/faI),r*(1-I(i-1)/c+ui)*(1.2-1.1*((mg*Ig+me*(1-Ig)))*I(i-1)/((n*(Q(i)-0.21*faqout*Q(i))*10000+0.5*wp*10000+0.8*P(i-1)*I(i-1)*2/3))));
    c = c*(1-0.003);%%土地盐碱率+Pop(i-1)
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
    if i+1957>2005
        ui = 0.01;
    end
    ki(i) = dI/I(i-1);
    %*(0.001*Pop(i-1)-0.01*UI)  *(0.001*Pop(i)) 0.0013*
    %*1.4*po/(Ig*I(i-1))
    I(i) = I(i-1)+dI;
     %% 水量
    if i+1957>1980
        elsewater = elsewater+0.006833;
    end
    costq(1) = 9.8250*10000*0.4;
    costq(i) = (1.1*mg*Ig+me*(1-Ig))*I(i);
    ds = (1.1*(mg*Ig+me*(1-Ig)))*I(i)/n-1*0.5*wp*10000/n-0.8*P(i-1)*I(i-1)*2/3/n;%这里0.8是纯井观区的取水，非纯灌区优先地下水
    %wmax(i) = min(ds,wp*10000);
    S(1) = 8.825*10000;
%     nn = 0;
%     while (Q(i)-0.15*faqout*Q(i))*10000-ds<0
%         ds=ds-0.3*wp*0.5*10000/n;
%         nn = nn+0.3;
%         if nn == 0.3
%             break
%         end
%     end
    downwater(i) = wp*(1);
    S(i) = elsewater*ds;
    dss(1) = 4.5*10000;
    dss(i) = max((Q(i)+0.21*faqout*Q(i))*10000-S(i),0);
    %% 分配给下游的水
    %dqout = qout(i-1)*min(0.01*(1-qout(i-1)/0.3/Q(i)/faqout),0.01*(1-qout(i-1)/0.3/Q(i))*(1-I(i)/c)*(1-1*(mg*Ig+me*(1-Ig))*I(i)/((n*(Q(i)-qout(i-1))*10000+0.6*wp*10000+0.6*2/30000*P(i)*I(i)))));
    %qout(i) = 0.15*faqout*Q(i);  
    %exps=[dI,costq(i)-costq(i-1),S(i)-S(i-1),downwater(i)-downwater(i-1)];
    %xlswrite('name.xlsx',exps,'Sheet2',num2str(i));
    end
    i = i+1;
end
%plot(years,I)
