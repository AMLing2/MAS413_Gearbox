% i_s1_min : minimum first stage gear ratio from Bild 15-38 lec 2 slide 12
% i_s1_max : maximum first stage gear ratio from Bild 15-38 lec 2 slide 12
% inc : step to increment between i_s1_min and i_s1_max
% z_1 : number of teeth in first stage pinion
% i_tot : total gear ratio requirement

% output format: [z_1,z_2,z_3,z_4,i_tot]
function [z_smallest,z_closest] = grat2stage(i_s1_min,i_s1_max,inc,z_1,i_tot)
  z_smallest = inf(1,5);
  z_closest = zeros(1,5);
  z_closest(5) = inf;
  s2_z1_min = 18; % same min as z_1 min
  s2_z1_max = 30; % arbitrary max, expensive to go over this, both min and max subject to change
  tol = 0.001; % only used for printing
  for s1_grat = i_s1_min:inc:i_s1_max
    z_2 = round(s1_grat * z_1);
    while gcd(z_2,z_1) > 1 % -1 to z2 until its relative prime (-1 gave closer final results)
      z_2 = z_2 - 1;
    end
    s2_grat = i_tot/s1_grat;
    for z_3 = s2_z1_min:s2_z1_max
      z_4 = round(s2_grat * z_3);
      while gcd(z_4,z_3) > 1 % same as above, +1 gave closer results here
        z_4 = z_4 + 1;
      end
      i_tot_t = (z_2/z_1) * (z_4/z_3); % calc total gear ratio

      % compare to smallest total z found:
      if ( abs(i_tot_t - i_tot)/i_tot ) < 0.01 % check if i_tot within 1% of i_tot
        if sum([z_1,z_2,z_3,z_4]) < sum(z_smallest(1:4)) % check if smaller total z than prev
            z_smallest(1) = z_1; z_smallest(2) = z_2; % save to array
            z_smallest(3) = z_3; z_smallest(4) = z_4;
            z_smallest(5) = i_tot_t;
        end
      end

      % compare to closest to i_tot found:
      if (abs(i_tot_t - i_tot)) < (abs(z_closest(5) - i_tot)) 
        z_closest(1) = z_1; z_closest(2) = z_2; % save to array
        z_closest(3) = z_3;  z_closest(4) = z_4;
        z_closest(5) = i_tot_t;
        if ismembertol(z_closest(5),i_tot,tol) % i_tot_o == i_tot
          sprintf("close match found: \n z_1: %d , z_2: %d \n z_3: %d , z_4: %d \n i_tot: %.6f ",z_1,z_2, z_3, z_4, i_tot_t)
        end
      end
    end
  end
 end
