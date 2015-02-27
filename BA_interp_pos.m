function g_ba_ci = BA_interp_pos( path )

addpath('./robotics3D/');

st = csvread([path, 'xkk_msckf.txt']);
states = reshape(st, [16, length(st)/16]);
g_p_i = states(14:16, :);
g_ba_ci = g_p_i(:,1);

for i = 2:length(g_p_i)
    g_p_cim5 = g_p_i(:,i-1);
    g_p_ci = g_p_i(:,i);
    
    g_p_cim4 = (g_p_cim5 * 4 + g_p_ci * 1) / 5;
    g_p_cim3 = (g_p_cim5 * 3 + g_p_ci * 2) / 5;
    g_p_cim2 = (g_p_cim5 * 2 + g_p_ci * 3) / 5;
    g_p_cim1 = (g_p_cim5 * 1 + g_p_ci * 4) / 5;
    
    g_ba_ci = [g_ba_ci, g_p_cim4, g_p_cim3, g_p_cim2, g_p_cim1, g_p_ci];
end

end
