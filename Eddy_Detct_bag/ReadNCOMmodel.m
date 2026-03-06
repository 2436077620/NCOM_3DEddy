function Model = ReadNCOMmodel(model_nm,res,hpkm)
    % 用于读取NCOM模式数据（主函数SWOTEddy_NCOMmodel.m）
%     res = 3.3; % 分辨率
    Model.lon = ncread(model_nm,'lon');
    Model.lat = ncread(model_nm,'lat');
    Model.depth = ncread(model_nm,'depth'); 
    Model.time = ncread(model_nm,'time')/24+datenum(2000,1,1); % 换算成以d为单位的数值
    Model.surf_el = ncread(model_nm,'surf_el');  % 海面高
    Model.water_u = ncread(model_nm,'water_u');  % 向东流速
    Model.water_v = ncread(model_nm,'water_v');  % 向北流速
    Model.salinity = ncread(model_nm,'salinity');% 盐度
    Model.temp = ncread(model_nm,'water_temp');  % 温度
    Model.water_w = ncread(model_nm,'water_w');  % 垂向流速
    Model.surf_atm_press = ncread(model_nm,'surf_atm_press');              % 海面大气压
    Model.surf_wnd_stress_gridx = ncread(model_nm,'surf_wnd_stress_gridx');% 海面向东风应力
    Model.surf_wnd_stress_gridy = ncread(model_nm,'surf_wnd_stress_gridy');% 海面向北风应力
    Model.surf_roughness = ncread(model_nm,'surf_roughness');              % 海面粗糙度
    Model.surf_temp_flux = ncread(model_nm,'surf_temp_flux');              % 海面温度谱
    Model.surf_solar_flux = ncread(model_nm,'surf_solar_flux');            % 海面短波谱
% ==========================温盐深：计算位势密度==============================
    SA = Model.salinity; % Absolute Salinity [g/kg] = [psu]
    CT = Model.temp;     % Conservative Temperature [deg C]
    [Y,~,Z] = meshgrid(Model.lat,Model.lon,Model.depth);
    p = ZP(Z,Y,0);       % sea pressure [dbar] 水深转为水压
    Model.potential_density = gsw_rho(SA,CT,p);
% ================对模式数据高通滤波，保留小于50km的成分======================
    Model.surf_el_hp = filt2(Model.surf_el,res,hpkm,'hp');     %分辨率res(3.3或8.3km)，滤波大小50km
    for i=1:length(Model.depth)                              %分层滤波
        Model.water_u_hp(:,:,i)  = filt2(Model.water_u(:,:,i),res,hpkm,'hp');
        Model.water_v_hp(:,:,i)  = filt2(Model.water_v(:,:,i),res,hpkm,'hp');
        Model.salinity_hp(:,:,i) = filt2(Model.salinity(:,:,i),res,hpkm,'hp');
        Model.temp_hp(:,:,i)     = filt2(Model.temp(:,:,i),res,hpkm,'hp');
        Model.potential_density_hp(:,:,i)= filt2(Model.potential_density(:,:,i),res,hpkm,'hp');
    end
end

