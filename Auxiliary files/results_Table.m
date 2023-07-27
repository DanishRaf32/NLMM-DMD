function results_Table(Error,FOMSnaptime,POD_ROMtime,NLMM_ROMtime,PODbasistime,NLMM_time,yNL,n_oc,rdefl,rdmd)

method={['FOM(n=',num2str(2*n_oc),')'];['POD-DEIM(r=',num2str(rdefl),')'];['NLMM (r=',num2str(rdefl),')',', DEIM(m=',num2str(rdmd),')']};
Offline_time={'-';FOMSnaptime;NLMM_time};
Online_time={'-';POD_ROMtime;NLMM_ROMtime};
Total_CPU_time={FOMSnaptime;POD_ROMtime+PODbasistime;NLMM_ROMtime+NLMM_time};
L1errorout1={'-';Error.err1PODout1Abs;Error.err1NLMMout1Abs};
L2errorout1={'-';Error.err2PODout1Abs;Error.err2NLMMout1Abs};
Linferrorout1={'-';Error.errInftyPODout1Abs;Error.errInftyNLMMOut1Abs};
Speed_ups={'-';FOMSnaptime/POD_ROMtime;FOMSnaptime/NLMM_ROMtime};

if size(yNL,2) >=2
    L1errorout2={'-';Error.err1PODout2Abs;Error.err1NLMMout2Abs};
    L2errorout2={'-';Error.err2PODout2Abs;Error.err2NLMMout2Abs};
    Linferrorout2={'-';Error.errInftyPODout2Abs;Error.errInftyNLMMOut2Abs};
    T=table(method,Offline_time,Online_time,Total_CPU_time,Speed_ups,L1errorout1,L2errorout1,Linferrorout1,L1errorout2,L2errorout2,Linferrorout2)
else
    T=table(method,Offline_time,Online_time,Total_CPU_time,Speed_ups,L1errorout1,L2errorout1,Linferrorout1)
end

end