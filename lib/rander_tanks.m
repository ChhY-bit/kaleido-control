function filepath = rander_tanks(x_data, u_data, params, dt, fps, filepath)
%RANDER_TANKS 根据四水箱数据和几何参数生成 GIF，不显示 figure 窗口。
%   filepath = rander_tanks(x_data,u_data,params,dt,fps,filepath)
%   x_data 为 4xN 水位 [m]，时刻为 (0:N-1)*dt [s]。
%   u_data 为 2xN 或 2x(N-1) 泵流量 [m^3/s]；后一种在末端保持最后输入。
%   params 使用 rpm.utils.load_quadruple_tank_params 返回的字段名。
%   fps 默认为 30；filepath 为文件夹或 .gif 文件名，默认 pwd/tanks.gif。
%   已有同名 GIF 在生成成功后覆盖；失败时保留原文件。
%   水位线性插值至帧时刻，输入零阶保持；末帧额外停留约 1/fps 秒。
%   GIF 延时按 0.01 s 量化，播放器实际播放速度可能有所不同。
%   水位须在箱高范围内；非负输入按最大流量限幅后映射总阀门转角。
%   gamma 在整段动画内保持不变，支管按初始 gamma 和相同设计流速定径。
%   仅需基础 MATLAB；函数通过 fprintf 输出进度，返回 GIF 完整路径。
%
%   rander_tanks(x_ext,u_ext,params,dt);
%   rander_tanks(x_ext,u_ext,params,dt,30,"tanks.gif");

    if nargin < 5 || isempty(fps), fps = 30; end
    if nargin < 6 || isempty(filepath), filepath = string(pwd); end
    validateattributes(dt, {'numeric'}, {'real','finite','scalar','positive'}, mfilename,'dt');
    validateattributes(fps, {'numeric'}, {'real','finite','scalar','positive','<=',100}, mfilename,'fps');
    validateattributes(x_data, {'numeric'}, {'real','finite','2d','nonempty','nrows',4,'nonnegative'}, mfilename,'x_data');
    validateattributes(u_data, {'numeric'}, {'real','finite','2d','nrows',2,'nonnegative'}, mfilename,'u_data');
    dt = double(dt); fps = double(fps);
    [area, hole, height, max_flow, speed, gamma] = read_geometry(params);
    n = size(x_data,2);
    if ~any(size(u_data,2) == [n, n-1]) || (n > 1 && isempty(u_data))
        error('rander_tanks:InputSize','u_data 须为 2xN 或 2x(N-1)。');
    end
    if any(any(x_data > height))
        error('rander_tanks:LevelRange','x_data 各水位须在 [0, 对应箱高] 内。');
    end
    if any(hole > area)
        error('rander_tanks:HoleSize','出水孔截面积不能超过对应水箱截面积。');
    end
    filepath = output_path(filepath);
    [output_dir,~,~] = fileparts(filepath);
    temporary = string(tempname(output_dir)) + ".gif";
    temporary_cleanup = onCleanup(@() remove_temporary(temporary)); %#ok<NASGU>

    duration = (n-1)*dt;
    if ~isfinite(duration)
        error('rander_tanks:Duration','数据总时长超出有限数值范围。');
    end
    frame_times = (0:floor(duration*fps))/fps;
    if duration > frame_times(end)
        frame_times(end+1) = duration;
    end
    % 累积量化，避免 30 fps 被逐帧固定舍入成 33.3 fps。
    centiseconds = round(100*[frame_times, duration+1/fps]);
    delays = max(1,diff(centiseconds))/100;
    frame_count = numel(frame_times);
    fprintf('开始生成四水箱 GIF：%.4g s，%d 帧，目标 %.4g fps。\n',duration,frame_count,fps);
    fprintf('保存位置：%s\n',filepath);
    if any(any(u_data > max_flow))
        fprintf('提示：超出泵最大流量的输入，仅在阀门显示时限幅。\n');
    end

    fig = figure('Visible','off','Color','w','Units','pixels', ...
        'Position',[0 0 1000 850],'MenuBar','none','ToolBar','none', ...
        'NumberTitle','off','Name','Four tanks','HandleVisibility','off');
    figure_cleanup = onCleanup(@() delete_figure(fig)); %#ok<NASGU>
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 10 8.5], ...
        'PaperSize',[10 8.5],'InvertHardcopy','off');
    ax = axes('Parent',fig,'Position',[.04 .04 .92 .88]);
    hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
    graphics = create_scene(ax,area,hole,height,max_flow,speed,gamma);
    heading = annotation(fig,'textbox',[.04 .93 .92 .05], ...
        'String','','EdgeColor','none','HorizontalAlignment','center', ...
        'FontSize',15,'FontWeight','bold','Interpreter','none');

    started = tic; last_report = -Inf;
    for k = 1:frame_count
        t = frame_times(k);
        position = min(n,1+t/dt);
        left = min(n,floor(position)); right = min(n,left+1);
        weight = position-left;
        levels = (1-weight)*double(x_data(:,left)) + weight*double(x_data(:,right));
        if isempty(u_data)
            flow = zeros(2,1); % 单个初始状态、没有输入区间
        else
            input_index = min(size(u_data,2),max(1,floor(t/dt+1e-9)+1));
            flow = min(double(u_data(:,input_index)),max_flow);
        end
        for i = 1:4
            bottom = graphics.bottom(i); w = graphics.width(i); c = graphics.center(i);
            set(graphics.water(i),'YData',[bottom bottom bottom+levels(i) bottom+levels(i)]);
            set(graphics.level(i),'String',sprintf('h%d = %.3f m',i,levels(i)));
            % 水体宽度固定为水箱直径，更新仅改变高度。
            set(graphics.water(i),'XData',c+[-w/2 w/2 w/2 -w/2]);
        end
        for j = 1:2
            angle = (pi/2)*flow(j)/max_flow(j);
            set_handle(graphics.main_handle(j),graphics.main_center(j,:),graphics.radius,angle);
        end
        set(heading,'String',sprintf('Quadruple tank  |  t = %.3f / %.3f s',t,duration));
        % print 直接获取隐藏 figure 的 RGB，不使用需要屏幕窗口的 getframe。
        rgb = print(fig,'-RGBImage','-r100');
        [indexed,map] = rgb2ind(rgb,256,'nodither');
        if k == 1
            imwrite(indexed,map,temporary,'gif','LoopCount',Inf,'DelayTime',delays(k));
        else
            imwrite(indexed,map,temporary,'gif','WriteMode','append','DelayTime',delays(k));
        end
        elapsed = toc(started);
        if k == 1 || k == frame_count || elapsed-last_report >= 2
            fprintf('进度 %6.2f%%  (%d/%d 帧)，已用 %.1f s，预计剩余 %.1f s\n', ...
                100*k/frame_count,k,frame_count,elapsed,elapsed/k*(frame_count-k));
            last_report = elapsed;
        end
    end
    [success,message] = movefile(temporary,filepath,'f');
    if ~success
        error('rander_tanks:SaveFailed','无法保存 GIF：%s',message);
    end
    fprintf('GIF 已保存：%s（耗时 %.1f s）。\n',filepath,toc(started));
end

function [area,hole,height,max_flow,speed,gamma] = read_geometry(params)
    validateattributes(params,{'struct'},{'scalar'},mfilename,'params');
    groups = {'TankArea',4;'HoleArea',4;'TankHeight',4; ...
              'PumpMaxFlow',2;'PumpFlowSpeed',2;'gamma',2};
    values = cell(6,1);
    for k = 1:6
        values{k} = zeros(groups{k,2},1);
        for i = 1:groups{k,2}
            name = sprintf('%s_%d',groups{k,1},i);
            if ~isfield(params,name)
                error('rander_tanks:MissingParameter','params 缺少字段 %s。',name);
            end
            if k == 6
                bounds = {'real','finite','scalar','>=',0,'<=',1};
            else
                bounds = {'real','finite','scalar','positive'};
            end
            validateattributes(params.(name),{'numeric'},bounds,mfilename,['params.' name]);
            values{k}(i) = double(params.(name));
        end
    end
    [area,hole,height,max_flow,speed,gamma] = values{:};
end

function filepath = output_path(filepath)
    if ~(ischar(filepath) && isrow(filepath)) && ~(isstring(filepath) && isscalar(filepath) && ~ismissing(filepath))
        error('rander_tanks:FilePath','filepath 须为字符行向量或字符串标量。');
    end
    filepath = string(filepath);
    if strlength(filepath) == 0, filepath = string(pwd); end
    [parent,name,extension] = fileparts(filepath);
    if isfolder(filepath) || (extension == "" && ~isfile(filepath))
        parent = filepath; name = "tanks"; extension = ".gif";
    elseif ~strcmpi(extension,'.gif')
        error('rander_tanks:FileExtension','文件名扩展名须为 .gif。');
    end
    if strlength(parent) == 0, parent = string(pwd); end
    if ~isfolder(parent)
        [success,message] = mkdir(parent);
        if ~success, error('rander_tanks:CreateDirectory','无法创建输出文件夹：%s',message); end
    end
    [~,attributes] = fileattrib(parent);
    filepath = string(fullfile(attributes.Name,name+extension));
end

function scene = create_scene(ax,area,hole,height,max_flow,speed,gamma)
    width = sqrt(4*area/pi); outlet = sqrt(4*hole/pi);
    main_width = sqrt(4*max_flow./speed/pi);
    scale = max([height;width;main_width]);
    gap = .85*scale;
    centers = [0; max(width([1 3]))/2+max(width([2 4]))/2+1.5*scale];
    center = centers([1 2 1 2]);
    bottom = [.65*scale;.65*scale; .65*scale+max(height(1:2))+gap; .65*scale+max(height(1:2))+gap];
    top = bottom+height;
    pump_x = [centers(1)-max(width([1 3]))/2-scale; centers(2)+max(width([2 4]))/2+scale];
    split_y = max(top(1:2))+.4*gap;
    radius = max(.13*scale,.65*max(main_width));
    % 回水池为固定示意布局，所有管径和水箱尺寸采用同一坐标比例。
    pool_left = pump_x(1)-.4*scale; pool_right = pump_x(2)+.4*scale;
    patch(ax,[pool_left pool_right pool_right pool_left],[0 0 .18 .18]*scale, ...
        [.80 .91 .98],'EdgeColor','none');
    plot(ax,[pool_left pool_left pool_right pool_right],[.3 0 0 .3]*scale,'Color',[.2 .3 .4],'LineWidth',1.4);
    for j = 1:2
        direction = 3-2*j; x = pump_x(j);
        pipe(ax,[x .16*scale; x split_y],main_width(j),scale);
        % 下层支路：分配 gamma，支管流速沿用总管设计流速。
        inlet_x = centers(j)-direction*.23*width(j);
        pipe(ax,[x split_y; inlet_x split_y; inlet_x top(j)+.06*scale], ...
            main_width(j)*sqrt(gamma(j)),scale);
        upper_tank = 5-j; % 泵 1 -> 水箱 4，泵 2 -> 水箱 3
        overhead = max(top(3:4))+(.65+.45*j)*scale;
        upper_inlet = center(upper_tank)+direction*.2*width(upper_tank);
        pipe(ax,[x split_y; x overhead; upper_inlet overhead; upper_inlet top(upper_tank)+.06*scale], ...
            main_width(j)*sqrt(1-gamma(j)),scale);
        pump_y = .45*scale;
        rectangle(ax,'Position',[x-radius pump_y-radius 2*radius 2*radius], ...
            'Curvature',[1 1],'FaceColor','w','EdgeColor',[.2 .3 .4]);
        patch(ax,x+radius*[-.55 .55 0],pump_y+radius*[-.45 -.45 .65], ...
            [.3 .5 .7],'EdgeColor','none');
        text(ax,x-direction*.32*scale,pump_y,sprintf('Pump %d',j), ...
            'HorizontalAlignment',alignment(j),'FontSize',10);
    end
    % 上层水箱排入对应下层水箱，下层排回水池。
    for i = 1:4
        if i <= 2, destination = .18*scale; else, destination = top(i-2)+.06*scale; end
        pipe(ax,[center(i) bottom(i); center(i) destination],outlet(i),scale);
    end
    scene.water = gobjects(4,1); scene.level = gobjects(4,1);
    for i = 1:4
        x = center(i); y = bottom(i); w = width(i); h = height(i);
        patch(ax,x+[-w/2 w/2 w/2 -w/2],[y y y+h y+h],'w','EdgeColor','none');
        scene.water(i) = patch(ax,x+[-w/2 w/2 w/2 -w/2],[y y y y], ...
            [.36 .73 .92],'EdgeColor','none');
        % 开口容器轮廓，底部按出水孔直径留出开口。
        plot(ax,x+[-w/2 -w/2 -outlet(i)/2],[y+h y y],'Color',[.15 .25 .35],'LineWidth',1.6);
        plot(ax,x+[outlet(i)/2 w/2 w/2],[y y y+h],'Color',[.15 .25 .35],'LineWidth',1.6);
        text(ax,x,y+h+.12*scale,sprintf('Tank %d',i),'HorizontalAlignment','center','FontWeight','bold');
        scene.level(i) = text(ax,x,y+.48*h,'','HorizontalAlignment','center','FontSize',10,'Color',[.1 .2 .3]);
    end
    scene.main_handle = gobjects(2,1); scene.main_center = zeros(2,2);
    for j = 1:2
        x = pump_x(j); y = .65*scale+.42*height(j);
        scene.main_center(j,:) = [x y];
        scene.main_handle(j) = valve(ax,[x y],radius,0,[.15 .4 .65]);
        text(ax,x+(3-2*j)*.22*scale,y-.28*scale,sprintf('u%d valve',j),'HorizontalAlignment','center','FontSize',9);
        valve(ax,[x split_y],radius,pi/2+(gamma(j)-.5)*pi/2,[.85 .48 .12]);
        text(ax,x,split_y+.25*scale,sprintf('gamma%d = %.2f',j,gamma(j)), ...
            'HorizontalAlignment','center','FontSize',9);
    end
    xlim(ax,[pool_left-.25*scale pool_right+.25*scale]);
    ylim(ax,[-.15*scale max(top(3:4))+2.1*scale]);
    scene.bottom = bottom; scene.center = center; scene.width = width; scene.radius = radius;
end

function value = alignment(j)
    if j == 1, value = 'right'; else, value = 'left'; end
end

function pipe(ax,points,width,scale)
%PIPE 通过数据坐标中的多边形表示真实管径，不用像素线宽代替几何宽度。
    color = [.70 .77 .82];
    if width == 0
        plot(ax,points(:,1),points(:,2),'--','Color',color,'LineWidth',.6);
        return;
    end
    for i = 1:size(points,1)-1
        a = points(i,:); b = points(i+1,:); v = b-a;
        if norm(v) <= eps(scale), continue; end
        normal = width/2*[-v(2) v(1)]/norm(v);
        corners = [a+normal;b+normal;b-normal;a-normal];
        patch(ax,corners(:,1),corners(:,2),color,'EdgeColor','none');
    end
    for i = 2:size(points,1)-1
        rectangle(ax,'Position',[points(i,:)-width/2 width width], ...
            'Curvature',[1 1],'FaceColor',color,'EdgeColor','none');
    end
end

function handle = valve(ax,center,radius,angle,color)
    rectangle(ax,'Position',[center-radius 2*radius 2*radius], ...
        'Curvature',[1 1],'FaceColor','w','EdgeColor',color,'LineWidth',1.2);
    handle = plot(ax,[0 0],[0 0],'-o','Color',color,'LineWidth',2,'MarkerSize',3);
    set_handle(handle,center,radius,angle);
end

function set_handle(handle,center,radius,angle)
    offset = .9*radius*[cos(angle) sin(angle)];
    set(handle,'XData',center(1)+[-offset(1) offset(1)], ...
        'YData',center(2)+[-offset(2) offset(2)]);
end

function delete_figure(fig)
    if isgraphics(fig), delete(fig); end
end

function remove_temporary(filepath)
    if isfile(filepath), delete(filepath); end
end
