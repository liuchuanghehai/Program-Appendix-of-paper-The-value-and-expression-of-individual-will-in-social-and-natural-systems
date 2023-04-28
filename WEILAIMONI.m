clc
clear
Q26=[8.7731
8.30298
10.34176
8.70939
8.13498
7.92453
7.41705
7.72178
9.54669
7.90031
8.32963
7.82972
9.72385
8.72433
7.43722
8.83091
10.1978
8.30666
9.46449
9.89454
8.82646
10.12786
9.89671
8.03253
8.34035
8.13338
8.59121
8.80016
9.43639
9.22616
8.63128
7.96613
7.34018
8.27427
8.52853
8.56756
9.39836
8.03581
8.60282
8.54215
7.79767
8.58083
8.98515
8.95105
8.34395
8.09371
8.6285
8.67111
8.36808
8.48605
8.91928
10.04437
8.7097
8.30514
9.00374
8.66779
7.85939
7.54524
9.622
8.42328
8.40409];
Q=[9.88245
8.69166
8.5067
9.12483
8.67744
8.31433
8.33036
8.22873
7.81298
7.86414
8.48668
8.62101
8.66644
8.6281
8.18904
8.46894
9.13813
8.4247
7.87377
8.83851
8.93991
9.07181
8.75247
9.33732
8.95495
8.11033
9.12703
10.74937
8.84667
9.57767
9.73591
10.84982
10.84071
9.64844
9.06194
8.79353
9.88762
9.50828
8.75362
8.4466
8.85022
9.7802
11.30721
9.14111
9.79077
9.44141
9.74536
11.2957
8.95507
9.74771
10.14172
8.90699
9.06026
8.98366
9.24925
8.23238
8.84699
8.9159
10.18195
9.17009
8.72015];
%wp_u = 4.18;
%wp_d = 0.86;
faI_m=0.6588;
faI_d=0.6588;
r = 0.1;
n = 0.6;
mg_m =350;
mg_d =480;
me_m =250;
me_d =250;
pop_d(1) = 26;
pop_m(1) = 104;
dpop = 0.003;
I_m=[];
I_m(1) =180 ;
I_d=[];
I_d(1) =87.45 ;
ui = 0;
Ig_m = 0.6;
Ig_d = 0.6;
wp_m = 4.18;
wp_d = 0.86;
years=[];
years(1) = 2020;
%% 情景
faqout = 1;
CI=[1,0];
WMAX=[1,1.5];
Cpop=[1.7,5;1.7,2.5];
QB = [0.2,0.4];
fo  = 1;
ds_m = 7*10000;
dw_d = 0.86*10000;
k=0.2;
%% 运算
i = 2;
while i+2019<=2080%从2021年开始运算，2020年给出
    years(i) = i+2019;
    P = 98+randn(1)*39.58;
    %% 中游
    cc_m = pop_m(i-1)*1.5;%Cpop;
    c_m = 300;
    faqout=faqout*(1-0.03);
    if faqout<0
        faqout=0;
    end
%     if i+2019>2030
%         faqout= 1;
%     end
    if c_m*faI_m<cc_m
        au = 0.03;
    else
        au = -0.03;
    end
    %au = 0.03;%工业强于农业此项不存在
    faI_m = faI_m*(1+au);
    if au>0
        faI_d = faI_d*(1+0.2*au);
    else
        faI_d = faI_d*(1);
    end
    if faI_m>0.8
        faI_m=0.8;
    end
    if faI_d>0.8
        faI_d=0.8;
    end
    k = k*(1+0.03);
    if k <=0.2
        k = 0.2;
    elseif k>0.5
        k=0.5;
    end 
    %k = 0.2;%工业强于农业此项不存在
    %wmax = 4.18;%*WMAX;
    wmax = 6;
    if wp_m <wmax
        if 10000*(Q(i)-2.9*faqout-k*Q(i))-ds_m>0         
            wp_m = wp_m+0;
        else
            wp_m = wp_m + min(0.3*0.4,wp_m*(1-wp_m/wmax+((I_m(i-1)*1.1*((mg_m*Ig_m+me_m*(1-Ig_m)))/(n*Q(i)*10000+0.6*10000*wp_m+0.6*P*I_m(i-1)*2/3)))));
        end
        if wp_m >=wmax
            wp_m =wmax;
        end
    else
        wp_m = wmax;
    end    
    dI_m = I_m(i-1)*min(r*(1-I_m(i-1)/(fo*c_m*faI_m+(1-fo)*cc_m)),r*(1-I_m(i-1)/c_m)*(1.2-1.1*((mg_m*Ig_m+me_m*(1-Ig_m)))*I_m(i-1)/((n*(Q(i)-2.9*faqout-k*Q(i))*10000+0.6*wp_m*10000+0.6*P*I_m(i-1)*2/3))));
    I_m(i) = I_m(i-1)+dI_m;
    pop_m(i) = pop_m(i-1)*(1+dpop)+ui*pop_d(i-1);    
%     if Q(i)>8.5
%         op =1.1 ;
%     elseif and(Q(i)>8,Q(i)<8.5)
%         op =1 ;
%     elseif and(Q(i)>7.5,Q(i)<8) 
         op =1 ;
%     elseif Q(i)<7.5
%         op =0.9 ;
%     end
    costq_m(1) = 4.585*10000;
    costq_m(i) = op*(1.1*mg_m*Ig_m+me_m*(1-Ig_m))*I_m(i);
    ds_m = costq_m(i)/(1-k)/n-0.6*wp_m*10000/n-0.6*P*I_m(i)*2/3/n;
    S(1) = 7*10000;
    downwater_m(1) = 4.18;
    downwater_m(i) = wp_m*(1);
    S(i) = min(ds_m,10000*(Q(i)-2.9*faqout));
    dss(1) = 3.9*10000;
    dss(i) = max((Q(i))*10000-S(i)+0.1*Q(i)*10000+(k-0.1)*0.3*Q(i)*10000,0);
    %% 下游
    cc_d = pop_d(i-1)*2.5;%Cpop;
    c_d = 220;   
    wmax_d = 6;%*WMAX; 

    if wp_d <wmax_d
        if dss*n+0.8*wp_d*10000+0.8*P*I_d(i-1)*2/3-1.1*(1.1*mg_d*Ig_d+me_d*(1-Ig_d))*I_d(i-1)*(I_d(i-1)+r*(1-I_d(i-1)/(fo*c_d*faI_d+(1-fo)*cc_d)))>0         
            wp_d = wp_d+0;
        else
            wp_d = wp_d + min(0.3*0.4,wp_d*(1-wp_d/wmax_d+((I_d(i-1)*1.1*((mg_d*Ig_d+me_d*(1-Ig_d)))/(n*dss(i)*10000+0.6*10000*wp_d+0.6*P*I_d(i-1)*2/3)))));
        end
        if wp_d >=wmax_d
            wp_d =wmax_d;
        end
    else
        wp_d = wmax_d;
    end
    dI_d = I_d(i-1)*min(r*(1-I_d(i-1)/(fo*c_d*faI_d+(1-fo)*cc_d)),r*(1-I_d(i-1)/c_d)*(1.2-1.1*((mg_d*Ig_d+me_d*(1-Ig_d)))*I_d(i-1)/((n*dss(i)+0.6*wp_d*10000+0.6*P*I_d(i-1)*2/3))));
    I_d(i) = I_d(i-1)+dI_d;
    pop_d(i) = pop_d(i-1)*(1+dpop-ui);  
    if pop_d(i)>27
        ui = 0.005;
    elseif pop_d(i)<25 
        ui = 0;
    end
    costq_d(1) = 4.585*10000;
    costq_d(i) = (1.1*mg_d*Ig_d+me_d*(1-Ig_d))*I_d(i);
    dw_d = 1.1*costq_d(i)/0.6-n*dss(i)/0.6-0.6*P*I_d(i)*2/3/0.6;
    downwater_d(1) = 0.86;
    downwater_d(i) = max(0,min(wp_d,dw_d));
%% 
    haveS_m(i) = Q(i)*10000-costq_m(i)-0.1*Q(i)*10000+P*I_m(i)-dss(i);
    haveS_d(i) = dss(i)-costq_d(i)+P*I_d(i);
    i = i+1;
end
plot(years,I_m)
hold on
plot(years,I_d)
