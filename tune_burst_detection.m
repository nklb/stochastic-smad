import burstDetection.generateSample;

% generate new sample of artificial data to test the performance
[trajecotries,labels]=generateSample(300);
len = size(trajecotries,1);
timepoints = 1:len;

%% plot examples from synthetic paths
figure(3);clf;
time = 0:5:1440;
for k = 1:9
    subplot(3,3,k);
    
    plot(time,trajecotries(:,k),'b-');hold on
    % plot labeled hills
    for l = 1:max(labels(:,k))
        take_this = labels(:,k)==l;
        plot(time(take_this),trajecotries(take_this,k),'r-','LineWidth',2);
    end
    title(sprintf('%d bursts',max(labels(:,k))));
    if ~isempty(intersect([1,4,7],k))
        ylabel('nuc/cyt SMAD2 ratio');
    end
    if ~isempty(intersect([7,8,9],k))
        xlabel('time / min');
    end
    
end

    
%% evaluete performance (ROC) of burst detection 

% load set of parameters to evaluate ROC-curve
burst_detection_params = load('data/ROC_parameters.mat','roc_performance','roc_parameters');
test_parameters = burst_detection_params.roc_parameters;    

%% run through the parameters and evaluete performance 
roc_test = zeros(size(test_parameters,1),2);
for i = 1:5:size(test_parameters,1)
    param = exp(test_parameters(i,:));
    roc_test(i,:) = eval_performance(param,trajecotries,timepoints,labels);
    fprintf('%3d: tpr = %3.0f %% ; fpr = %3.0f %%\n',i,roc_test(i,:).*100);
end

figure;
plot(roc_test(:,2),roc_test(:,1),'ko');
xlabel('False positve rate');
ylabel('True positive rate');
title('Burst detection ROC curve')
hold on;

%% evaluate performance of default parameters
tpr_fpr=eval_performance([],trajecotries,timepoints,labels);
plot(tpr_fpr(2),tpr_fpr(1),'r*');


function ret = eval_performance(param,trajecotries,timepoints,label)
    import burstDetection.burstDetect;
    % if parameters are not defined, use default paramters
    if isempty(param)
        [select_as_burst,~]=burstDetect(trajecotries,timepoints);
    else
        [select_as_burst,~]=burstDetect(trajecotries,timepoints,log(param));
    end

    % True positive rate
    hit=select_as_burst>0 & label>0;
    TPR = mean(sum(hit)./sum(label>0));

    % False positive rate
    false_predict = select_as_burst>0 & label==0;
    FPR = mean(sum(false_predict)./sum(label==0));

    ret=[TPR,FPR];
end
