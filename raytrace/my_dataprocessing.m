function my_dataprocessing(path,proportional_scale,varagin)   

%% default value
condition_num = 21;
keyword = '后';

%% methods


%% main 

close all
save_bubble =zeros(condition_num,12);
for k = 1 : condition_num
    close all
    clearvars -except k save_bubble path condition_num keyword proportional_scale;
    ccc = (k-1)*100;
    file_path = [path '\' num2str(ccc) '\']; 
    if exist([path '\' num2str(ccc)],'dir') ~= 7
        warning([path '\' num2str(ccc) '文件夹不存在！']);
        continue
    end
    getfilename1 = ls([file_path,'*.xlsx']');             
    filename1 = cellstr ( getfilename1 ) ;   
    [L1,~]=size(getfilename1);
    indata = cell(L1,2);
    
    E_num = 0;

    for i=1:L1
        filename_cur = char(filename1(i));
        if ~isempty(strfind(filename_cur,keyword)) & isempty(strfind(filename_cur,'~'))
            E_num = E_num + 1;
            filename_full = [file_path, filename_cur];
            indata{E_num,1} = xlsread(filename_full);
            indata{E_num,2} = indata{E_num,1}(:,1); 
        end
    end
    
    
    bubble = zeros();
    bubble_num = 0;
    dia_3 = 0;
    dia_2 = 0;
    for i = 1 : E_num
        [aaa,~] = size(indata{i,2});
        for j = 1 : aaa
            bubble_num = bubble_num + 1;
            bubble(bubble_num,1) = indata{i,2}(j,1) * proportional_scale;
            dia_3 = dia_3 + bubble(bubble_num,1)^3;
            dia_2 = dia_2 + bubble(bubble_num,1)^2; 
        end
    end

    for i = 1 : bubble_num
        if bubble(i,1) < 0.3
            save_bubble(k,5) = save_bubble(k,5) + 1;
        elseif bubble(i,1) < 0.6
            save_bubble(k,6) = save_bubble(k,6) + 1;
        end
        if bubble(i,1) > 0.9
            save_bubble(k,7) = save_bubble(k,7) + 1;
        end
        if bubble(i,1) > 1.2
            save_bubble(k,8) = save_bubble(k,8) + 1;
        end
        if bubble(i,1) > 1.5
            save_bubble(k,9) = save_bubble(k,9) + 1;
        end
        if bubble(i,1) > 1.8
            save_bubble(k,10) = save_bubble(k,10) + 1;
        end
        if bubble(i,1) > 2
            save_bubble(k,11) = save_bubble(k,11) + 1;
        end
        
    end
    save_bubble(k,5) = save_bubble(k,5)/bubble_num;
    save_bubble(k,6) = save_bubble(k,6)/bubble_num;
    save_bubble(k,7) = save_bubble(k,7)/bubble_num;
    save_bubble(k,8) = save_bubble(k,8)/bubble_num;
    save_bubble(k,9) = save_bubble(k,9)/bubble_num;
    save_bubble(k,10) = save_bubble(k,10)/bubble_num;
    save_bubble(k,11) = save_bubble(k,11)/bubble_num;


    save_bubble(k,1) = mean(bubble);
    save_bubble(k,2) = dia_3/dia_2;
    save_bubble(k,3) = max(bubble);
    save_bubble(k,4) = min(bubble);
    
    save_bubble(k,12) = (bubble_num/E_num)*100;

    %作图
    figure
    h2 = histogram(bubble,'Normalization','probability','BinWidth',0.05);
    if exist([path '\result'],'dir') ~= 7
        mkdir([path '\result'])
    end
    saveas(gcf,[path '\result\',num2str(ccc),'.bmp']);
end
xlswrite([path '\result\result.xlsx'],save_bubble);
    
end