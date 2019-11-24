function R = rot_world2local(ini_normal)

	Dp = dot(ini_normal, [0;0;1]);
    Cp = cross(ini_normal, [0;0;1]);
    I = eye(3);
    Cp_x = [0 -Cp(3) Cp(2); Cp(3) 0 -Cp(1); -Cp(2) Cp(1) 0];
    R = (Cp_x*Cp_x)/(1+Dp) + Cp_x + I;

end