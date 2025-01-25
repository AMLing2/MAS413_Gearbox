function [z_1_o,z_2_o,z_3_o,z_4_o,i_tot_o] = grat2stage(i_s1_min,i_s1_max,inc,z_1,i_tot)
  z_1_o = 0; z_2_o = 0; z_3_o = 0; z_4_o = 0;
  i_tot_o = inf;
  s2_z1_min = 18; % same min as z_1 min
  s2_z1_max = 30; % arbitrary max, probably bad to go over this, both min and max subject to change
  tol = 0.001;
  for s1_grat = i_s1_min:inc:i_s1_max
    z_2 = round(s1_grat * z_1);
    while gcd(z_2,z_1) > 1 % -1 to z2 until its relative thingy (-1 gave better results)
      z_2 = z_2 - 1;
    end
    s2_grat = i_tot/s1_grat;
    for z_3 = s2_z1_min:s2_z1_max
      z_4 = round(s2_grat * z_3);
      while gcd(z_4,z_3) > 1 % same as above, +1 gave better results here...
        z_4 = z_4 + 1;
      end
      i_tot_t = (z_2/z_1) * (z_4/z_3); % calc total gear ratio
      if (abs(i_tot_t - i_tot)) < (abs(i_tot_o - i_tot)) % compare to lowest found
        i_tot_o = i_tot_t;
        z_1_o = z_1; z_2_o = z_2;  z_3_o = z_3;  z_4_o = z_4;
        if ismembertol(i_tot_o,i_tot,tol) % i_tot_o == i_tot
          sprintf("close match found: \n z_1: %d , z_2: %d \n z_3: %d , z_4: %d \n i_tot: %.6f ",z_1,z_2, z_3, z_4, i_tot_o)
        end
      end
    end
  end
 end
