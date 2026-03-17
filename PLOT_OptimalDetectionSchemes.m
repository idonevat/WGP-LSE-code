
clear; clc;
close all;
rng(11);


warps = { ...
    make_warp_beta(0.5, 0.5), ...
    make_warp_gamma(0.5, 1.0), ...
    make_warp_weibull(1.0, 5.0) ...
    };

warpNames = { 'Beta ', 'Gamma ', 'Weibull ' };

marker_shapes = {'o', '+', '*', 's', 'd'};

%% -------------------------
% Run per warp
%% -------------------------
figure('Name','Sensitivity','Color','w','Position',[100 100 1200 850])
t = tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
t.OuterPosition = [0 0.08 1 0.92];   % leave space for legend

% handles for one global legend
hLegend = gobjects(5,1);

for iw = 1:numel(warps)
    name_WS=['WS_optimal_detection_schemes_' lower(warpNames{iw}),'.mat'];
    load(name_WS)

    %% -------------------------
    % Plot: 2-panel per warp
    %% -------------------------
    % (Left) J1 for R1 (moments objective)
    nexttile;
    for ie = 1:numel(ellList)
        y_smooth = smooth(Nlist, J1_R1(ie,:), 0.2, 'loess');
        current_marker = marker_shapes{mod(ie-1, length(marker_shapes)) + 1};
        %h= plot(Nlist, J1_R1(ie,:), 'Marker', current_marker, 'LineWidth', 2,'MarkerSize', 10);
        h= plot(Nlist, y_smooth, 'Marker', current_marker, 'LineWidth', 1,'MarkerSize', 6);
        h.MarkerIndices = 1:3:length(Nlist);
        hold on;
    end
    grid on;
    xlabel('Number of observations','FontWeight','bold','FontSize',16);
    %ylabel('Spatial average moments-risk','FontWeight','bold','FontSize',16);
    ylabel('Moments SIE','FontWeight','bold','FontSize',16);
    %title([warpNames{iw} 'process : moments-optimal rule'],'FontWeight','bold','FontSize',16);
    title([warpNames{iw} ' process: moments-based rule'],'FontWeight','bold','FontSize',16);
    %legend(arrayfun(@(e) sprintf('\\ell=%.2f', e), ellList, 'UniformOutput', false), ...
    %    'Location','northeast' ,'Interpreter','latex');
    %[hleg, hobj] = legend;
    xticks([Nlist(1):50:Nlist(end)])
    %set(hleg,'Location','best','NumColumns',1,'FontSize',20);
    xlim([Nlist(1) Nlist(end)])
    set(gca,'FontSize',12)
    set(gca,'Fontweight','bold')
    % (Right) J2 for R2 (exceedance objective)
    nexttile;
    for ie = 1:numel(ellList)
        y_smooth = smooth(Nlist, J2_R2(ie,:), 0.2, 'loess');
        current_marker = marker_shapes{mod(ie-1, length(marker_shapes)) + 1};
        h= plot(Nlist, y_smooth, 'Marker', current_marker, 'LineWidth', 1,'MarkerSize', 6);
        h.MarkerIndices = 1:3:length(Nlist);
        hold on;
        if iw == 1
            hLegend(ie) = h;

        end
    end
    grid on;
    xlabel('Number of observations','FontWeight','bold','FontSize',16);
    %ylabel('Spatial average exceedance-risk ','FontWeight','bold','FontSize',16);
    ylabel('Exceedance SIE','FontWeight','bold','FontSize',16);
    title([warpNames{iw} ' process: exceedance-based rule'],'FontWeight','bold','FontSize',16);
    %legend(arrayfun(@(e) sprintf('\\ell=%.2f', e), ellList, 'UniformOutput', false), ...
    %    'Location','northeast' ,'Interpreter','latex');
    %[hleg, hobj] = legend;
    xticks([Nlist(1):50:Nlist(end)])
    %set(hleg,'Location','best','NumColumns',1,'FontSize',20);
    xlim([Nlist(1) Nlist(end)])
    set(gca,'FontSize',12)
    set(gca,'Fontweight','bold')

    % if iw == 1
    %     hLegend(2*iN-1) = h;
    %     hLegend(2*iN)   = h;
    % end


end
%% ============================================================
% GLOBAL LEGEND
%% ============================================================
lgd = legend(hLegend, ...
    arrayfun(@(n) sprintf('l=%-10.2f',n), ellList, 'UniformOutput', false), ...
    'Orientation','horizontal', ...
    'FontSize',14, ...
    'NumColumns',5);

lgd.Units = 'normalized';
lgd.Position = [0.15 -0.08 0.6 0.07];
pause(1)
exportgraphics(gcf, ['perfect_match_all.png'], 'Resolution', 300);

