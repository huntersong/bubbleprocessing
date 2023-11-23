function my_histogram(Indata,varargin)
% 输入参数Indata为待画直方图的数据，格式为N×1的数组，默认只画'pdf'概率密度图
% x轴为直径（mm），y轴为pdf
% bincolor  'blue' [0.6902 0.7765 0.9020]   'green' [0.6745,0.8431,0.6706] 'red' [0.9765 0.7216 0.7137]
% linecolor  'blue' [0 0.4431 0.7412]  'green'[0.1137 0.6902 0.2706] 'red' [0.6863 0.0510 0.0824]

%% default value
data_max = max(Indata);
binwidth = round(data_max)/50;
edges = 0:binwidth:(fix(data_max)+1);
pHat = lognfit(Indata);
y = lognpdf(edges,pHat(1),pHat(2));

bincolor = [0.6902 0.7765 0.9020
         0.6745,0.8431,0.6706
         0.9765 0.7216 0.7137];
linecolor = [0.0000 0.4431 0.7412
             0.1137 0.6902 0.2706
             0.6863 0.0510 0.0824];
B = 2; L = 3;  % 直方图的默认颜色是绿色，线的默认颜色是红色，均为半透明
saveflag = 0;
fitflag = 1;
y_min = 0;
y_max = pdf_max(Indata,edges);

if numel(varargin) > 0
    tf = {'save','bincolor','linecolor','no_log_fit'};
    for i = 1 : numel(tf)
        for j = 1 : numel(varargin)
            tf_1 = strcmp(tf{i},varargin{j});
            if tf_1 == 1
                switch i
                    case 1
                        saveflag = 1;
                    case 2
                        B = varargin{j+1};
                    case 3
                        L = varargin{j+1};
                    case 4
                        fitflag = 0;
                end
            end
        end
    end
end


%% main

if saveflag == 0
    figure
    histogram(Indata,'BinEdges',edges,'Normalization','pdf','EdgeAlpha',0.5,'EdgeColor','k','FaceAlpha',0.5,...
        'FaceColor',bincolor(B,:),'LineWidth',0.5,'Parent',app.UIAxes_bubblesize);
    hold on
    if fitflag == 1
        plot1 = plot(edges,y,'Linewidth',1,'color',linecolor(L,:));
        plot1.Color(4) = 0.67;  %设置拟合曲线的透明度
    end
    axis([0 (fix(data_max)+1) y_min round(10^2*y_max)/10^2]);
    set(gca,'FontSize',17,'Fontname','Arial');
    legend('PDF of bubble size','Lognormal fit of data');
    yticks([0:0.1:round(10^2*y_max)/10^2]);		%具体的y轴刻度
    pos=axis;%取得当前坐标轴的范围，即[xmin xmax ymin ymax]
    xlabel('Bubble Diameter(mm)','fontsize',20,'fontname','Arial','position',[(pos(1)+pos(2))/2 pos(3)-0.06])	%设置x轴字体
    ylabel('Probability Density Function','fontsize',20,'fontname','Arial','position',[pos(1)-0.2 (pos(3)+pos(4))/2 ])	%设置y轴字体
    hold off
else
    figure
    histogram(Indata,'BinEdges',edges,'Normalization','pdf','EdgeAlpha',0.5,'EdgeColor','k','FaceAlpha',0.5,...
        'FaceColor',bincolor(B,:),'LineWidth',0.5);
    hold on
    if fitflag == 1
        plot1 = plot(edges,y,'Linewidth',1,'color',linecolor(L,:));
        plot1.Color(4) = 0.67;  %设置拟合曲线的透明度
    end
    axis([0 (fix(data_max)+1) y_min round(10^2*y_max)/10^2]);
    set(gcf,'position',[400,150,960,740]);%设置画图的大小
    set(gca,'FontSize',17,'Fontname','Arial');
    legend('PDF of bubble size','Lognormal fit of data');
    yticks([0:0.1:round(10^2*y_max)/10^2]);		%具体的y轴刻度
    pos=axis;%取得当前坐标轴的范围，即[xmin xmax ymin ymax]
    xlabel('Bubble Diameter(mm)','fontsize',20,'fontname','Arial','position',[(pos(1)+pos(2))/2 pos(3)-0.06])	%设置x轴字体
    ylabel('Probability Density Function','fontsize',20,'fontname','Arial','position',[pos(1)-0.2 (pos(3)+pos(4))/2 ])	%设置y轴字体
    savepath = uigetdir({},'选择文件夹');
    filename = inputdlg('输入图片名称：');
    print('-dpng',[savepath '\' filename{1,1} '.png']); %保存图像
    hold off
end

end


function y_max = pdf_max(Indata,edges)
N = histcounts(Indata,'BinEdges',edges,'Normalization','pdf');
y_max = max(N);
end