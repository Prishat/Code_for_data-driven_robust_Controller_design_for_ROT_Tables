function [C, Ta] = heatTranStep(q, T, Taexp, Ta0, Ta1, dt, varC, varTa, varlik, Nparticlesamp)

ns = size(T,1);
alpha = 1;  beta = 0;
Ta = [Ta0; Ta1; zeros(ns-2,1)];
C = zeros(ns,1); C(1) = (T(2)-T(1))./(Ta(1)-T(1))./dt(1);
deltaTa = Ta1-Ta0;
dq = q(2) - q(1);
% Ta(2) = Ta(1) + deltaTa;


% % Lognormal mean, var from Gaussian mean, var parameters
% meanln = @(mu,sig) exp(mu + sig);
% varln = @(mu,sig) exp(2.*mu + sig).*(exp(sig)-1);
% 
% priorC = @(cm1,c) pdf('Lognormal',c,meanln(cm1,varC),varln(cm1,varC));
% priorTa = @(Ta,Tap1) pdf('Lognormal',Tap1,meanln(Ta,varTa),varln(Ta,varTa));   % pdf('lognormal',Ta-Tap1+1,1,varTa);




% Gaussian mean, var from Lognormal mean, var parameters
meanG = @(mu,sig) log(mu./sqrt(1+sig./mu.^2));
varG = @(mu,sig) log(1+sig./mu.^2);

priorC = @(cm1,c)    pdf('Normal', log(c),    meanG(cm1,varC), varG(cm1,varC));
priorTa = @(Ta,Tap1) pdf('Normal', log(Tap1), meanG(Ta,varTa), varG(Ta,varTa));





modQTPlus1 = @(q,T,c,Ta)  q - c.*(Ta-T);            % to be evaluated at ii+1
liklihood = @(q,T,c,Ta) pdf('Normal',((modQTPlus1(q,T,c,Ta))),0,varlik);

% randgenC = @(cm1,varC,n) pdf('Normal',c,cm1,varC,n);
% randgenTa = @(Ta,varTa,n) pdf('Normal',Tap1,Ta,varTa,n);

sigmaCTafull = zeros(2,2,ns);
sigmaCTafull(:,:,1) = [varC,0;0,varTa];

% figure()
% plot(1:ns, Taexp, '.b'); hold on;


for ii=2:ns-2

    detC = (T(ii+1)-T(ii))./(Ta(ii)-T(ii))./dt(ii);

    detTaPlus1 = Ta(ii) + ...
        alpha.*deltaTa./(T(ii)-T(ii-1)).*(T(ii+1) - T(ii)) + ...
        beta.*deltaTa./(q(ii)-q(ii-1)).*(q(ii+1) - q(ii));


    [postPDFmean, postPDFvar] = jointpostPDF(liklihood, priorC,priorTa, ...
        q(ii+1),T(ii+1),[C(ii-1),detC],[Ta(ii),detTaPlus1], ...
        sigmaCTafull(:,:,ii-1), Nparticlesamp, modQTPlus1);

    % update prior variance with new posterior variance
    sigmaCTafull(:,:,ii) = [varC,0;0,varTa];  %postPDFvar;
    C(ii) = postPDFmean(1);
    Ta(ii+1) = postPDFmean(2);
    deltaTa = Ta(ii+1) - Ta(ii);


    % plot(ii, Ta(ii+1), '.r');
    % plot(ii, Ta(ii+1)+1e2*postPDFvar(2,2), '-k')
    % plot(ii, Ta(ii+1)-1e2*postPDFvar(2,2), '-k')


    if ii==9
       disp('The 10th step...')
    end

end


end







function [meanC_Ta, sigC_Ta] = jointpostPDF(liklihood, ...
    priorC, priorTa, q,T,c,Ta, Sigma_C_Ta, Nsamp, modQTPlus1)

% Gaussian mean, var from Lognormal mean, var parameters
meanG = @(mu,sig) log(mu./sqrt(1+sig./mu.^2));
varG = @(mu,sig) log(1+sig./mu.^2);

jointMu = [c(2) Ta(1)];
sig = Sigma_C_Ta;
% R = chol(sig);
% C_Ta = repmat(jointMu,Nsamp,1) + randn(Nsamp,2)*R;  % Nsamp X 2 rand matrix
try
    % C_Talog = lhsnorm(log(jointMu),sig,Nsamp,'off');
    C_log  = lhsnorm(meanG(jointMu(1),sig(1,1)),varG(jointMu(1),sig(1,1)),Nsamp,'off');
    Ta_log = lhsnorm(meanG(jointMu(2),sig(2,2)),varG(jointMu(2),sig(2,2)),Nsamp,'off');
    C_Talog = [C_log, Ta_log];
    C_Ta = exp(C_Talog);
    % C_Ta = mvnrnd(jointMu,sig,Nsamp);
catch
    error('An unexpected error has occured')
end


% figure()
% plot(C_Ta(:,1), C_Ta(:,2), '.b')

Cpt = sort(C_Ta(:,1));
Tapt = sort(C_Ta(:,2));
[Cmesh,Tamesh] = meshgrid(Cpt,Tapt);
Csort = Cmesh(:); Tasort = Tamesh(:);

dc = abs(Cpt(2:end) - Cpt(1:end-1));
dc = [dc; dc(1)];
dTa = abs(Tapt(2:end) - Tapt(1:end-1));
dTa = [dTa; dTa(1)];
[dcMesh,dTaMesh] = meshgrid(dc,dTa);

likval = liklihood(q,T,Csort,Tasort);

% figure()
% plot(Csort,pC./sum(pC),'-b',Csort,likval./sum(likval),'-r', ...
%     Csort, pC.*likval./sum(pC.*likval),'-k')

% figure()
% plot(Tasort,pTa./sum(pTa),'-b',Tasort,likval./sum(likval),'-r', ...
%     Tasort, pTa.*likval./sum(pTa.*likval),'-k')

% figure()
% scatter3(Csort, Tasort, likval,4,likval)

% pC = priorC(c(1),Csort).*priorC(c(2),Csort);
% pTa = priorTa(Ta(1),Tasort).*priorTa(Ta(2),Tasort);
pC = priorC(c(2),Csort);
pTa = priorTa(Ta(1),Tasort);

% figure()
% plot(linspace(c(1)*0.5,c(1)*1.5,1000), ...
%     priorC(c(1),linspace(c(1)*0.5,c(1)*1.5,1000))); hold on
% plot(Csort,priorC(c(1),Csort),'.r'); 
% plot(Csort,priorC(c(2),Csort),'.k'); 
% hold off
% 
% figure()
% plot(linspace(Ta(1)*0.5,Ta(1)*1.5,1000), ...
%     priorTa(Ta(1),linspace(Ta(1)*0.5,Ta(1)*1.5,1000))); hold on
% plot(Tasort,priorTa(Ta(1),Tasort),'.r'); 
% plot(Tasort,priorTa(Ta(2),Tasort),'.k'); 
% hold off


% the following equation can take combination of scalar and vectors
pdfsamp = likval.*pC.*pTa.*dcMesh(:).*dTaMesh(:);
pdfsampnorm = pdfsamp./sum(pdfsamp);



% 
% errmod = modQTPlus1(q,T,Csort,Tasort);
% pdferr = pdf('Normal',errmod,0,3e0);
% figure()
% scatter3(Csort,Tasort,likval,'.b')





% pdfsamp = liklihood(q,T,C_Ta(:,1),C_Ta(:,2));

meanC = sum(Csort.*pdfsampnorm);
meanTa = sum(Tasort.*pdfsampnorm);

varC = sum((Csort-meanC).^2.*pdfsampnorm);
varTa = sum((Tasort-meanTa).^2.*pdfsampnorm);
varC_Ta = sum((Tasort-meanTa).*(Csort-meanC).*pdfsampnorm);


[~,maxCInd] = max(Csort.*pdfsampnorm);
[~,maxTaInd] = max(Tasort.*pdfsampnorm);
% meanC_Ta = [Csort(maxCInd), Tasort(maxTaInd)];
meanC_Ta = [meanC, meanTa];
sigC_Ta = [varC, varC_Ta; varC_Ta, varTa];

end



% pCsort = pdf('Normal',Csort,c(1),1e-3);
% figure()
% plot(Csort,pCsort,'.b')


% figure()
% scatter3(Csort, Tasort, pdfsampnorm,4,pdfsampnorm)


% function [postPDFmean, postPDFvar, postPDF] = jointpostPDF(modQTPlus1, ...
%     liklihood,priorC,priorTa, q,T,c,Ta, varC, varTa, varlik, Ns)
% 
%     errQ = modQTPlus1(q,T,c(2),Ta(2));
%     postPDFvar = 1./(1./varC + 1./varTa + 1./varlik);
%     postPDFmean = postPDFvar.*(errQ.*0./varlik + c(1)./varC + Ta(1)./varTa);  % 0 with varlik since mean is zero
%     postPDF = liklihood(q,T,c(2),Ta(2)).*priorC(c(1),c(2)).*priorTa(Ta(1),Ta(2));
% 
% end




