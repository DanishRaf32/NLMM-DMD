function [Error,quality] = error_cal(yFOM,yROM_POD,yROM_NLMM)

%% prerequisites
Error.errorNLMM     = abs(yFOM-yROM_NLMM);
Error.errorNLMMOut1 = abs(Error.errorNLMM(:,1));
Error.errorNLMMOut2 = abs(Error.errorNLMM(:,2));

Error.errorPOD      = abs(yFOM-yROM_POD);
Error.errorPODout1  = abs(Error.errorPOD(:,1));
Error.errorPODout2  = abs(Error.errorPOD(:,2));

%% L1-norm
% absolute Error
Error.err1NLMMout1Abs = norm(Error.errorNLMMOut1,1); 
Error.err1NLMMout2Abs = norm(Error.errorNLMMOut2,1);

Error.err1PODout1Abs = norm(Error.errorPODout1,1);
Error.err1PODout2Abs = norm(Error.errorPODout2,1);

% relative Error
Error.err1NLMMout1Rel = norm(Error.errorNLMMOut1,1)./norm(yFOM(1,:),1); % err1 = sum(abs(error))
Error.err1NLMMout2Rel = norm(Error.errorNLMMOut2,1)./norm(yFOM(2,:),1);

Error.err1PODout1Rel = norm(Error.errorPODout1,1)./norm(yFOM(1,:),1);
Error.err1PODout2Rel = norm(Error.errorPODout2,1)./norm(yFOM(2,:),1);

%% L2-norm
% absolute Error
Error.err2NLMMout1Abs = norm(Error.errorNLMMOut1,2); 
Error.err2NLMMout2Abs = norm(Error.errorNLMMOut2,2);

Error.err2PODout1Abs = norm(Error.errorPODout1,2);
Error.err2PODout2Abs = norm(Error.errorPODout2,2);


% relative Error
Error.err2out1 = norm(Error.errorNLMMOut1,2)./norm(yFOM(1,:),2); % err2 = sqrt(error'*error)
Error.err2out2 = norm(Error.errorNLMMOut2,2)./norm(yFOM(2,:),2);

Error.err2PODout1 = norm(Error.errorPODout1,2)./norm(yFOM(1,:),2);
Error.err2PODout2 = norm(Error.errorPODout2,2)./norm(yFOM(2,:),2);

%% Linf-norm
% absolute Error
Error.errInftyNLMMOut1Abs = norm(Error.errorNLMMOut1,Inf);
Error.errInftyNLMMOut2Abs = norm(Error.errorNLMMOut2,Inf);

Error.errInftyPODout1Abs = norm(Error.errorPODout1,Inf);
Error.errInftyPODout2Abs = norm(Error.errorPODout2,Inf);


% relative Error 
Error.errInftyNLMMOut1 = norm(Error.errorNLMMOut1,Inf)./norm(yFOM(1,:),Inf); % errInfty = max(abs(error))
Error.errInftyNLMMOut2 = norm(Error.errorNLMMOut2,Inf)./norm(yFOM(2,:),Inf);

Error.errInftyPODout1 = norm(Error.errorPODout1,Inf)./norm(yFOM(1,:),Inf);
Error.errInftyPODout2 = norm(Error.errorPODout2,Inf)./norm(yFOM(2,:),Inf);

%% Quality measure 
% absolute 
quality.qm1NLMM = Error.err1NLMMout1Abs + Error.err1NLMMout2Abs;
quality.qm1POD = Error.err1PODout1Abs + Error.err1PODout2Abs;

quality.qm2NLMM = Error.err2NLMMout1Abs + Error.err2NLMMout2Abs;
quality.qm2POD = Error.err2PODout1Abs + Error.err2PODout2Abs;

quality.qmInftyNLMM = Error.errInftyNLMMOut1Abs + Error.errInftyNLMMOut2Abs;
quality.qmInftyPOD = Error.errInftyPODout1Abs + Error.errInftyPODout2Abs;
    
% relative 
% % qm1LPTD = err1out1 + err1out2
% % qm1POD = err1PODout1 + err1PODout2
% % 
% % qm2LPTD = err2out1 + err2out2
% % qm2POD = err2PODout1 + err2PODout2
% % 
% % qmInftyLPTD = errInftyOut1 + errInftyOut2
% % qmInftyPOD = errInftyPODout1 + errInftyPODout2
end