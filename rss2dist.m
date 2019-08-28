function [dist_RSS] = rss2dist(RSS,NLOS, d0, RU, freq)
%% Converting RSS to distance:
% ----------------------------------------------
% INPUTS:
% NAME:         SIZE:       DESCRIPTION:
% RSS           vector      
% NLOS          vector      RSS NLOS detection result.      
% d0            1x1         Reference distance in the path-loss model.
%                           (It is important that the same value is used 
%                           when estimating the path-loss parameters.)
% RU            struct      RU definitions.
% freq          str         RSS frequency currently used.
% ----------------------------------------------
% OUTPUT:
% NAME:         SIZE:       DESCRIPTION:
% dist_RSS      vector      RSS for all RUs converted to distance.

for nn = 1:size(NLOS,2)

    if NLOS(nn) == 0 % If NLOS is NOT detected; use LOS model:
        
        if strcmp(freq,'24')
            dist_RSS(1,nn) = d0*10.^((RSS(1,nn) - RU(nn).Pt - RU(nn).PL024_LOS)./(-10*RU(nn).n24_LOS));
        elseif strcmp(freq,'58')
            dist_RSS(1,nn) = d0*10.^((RSS(1,nn) - RU(nn).Pt - RU(nn).PL058_LOS)./(-10*RU(nn).n58_LOS));
        end
        
    else % If NLOS is detected;
        
        % Use LOS model from RU to reference (reflection/diffraction) point:
        % (This is the case because the MU is measuring.)
        dist_RU_refl = dist(RU(nn).position, RU(nn).NLOSrefl');
       
        % Use NLOS model from reference (reflection/diffraction) point to MU:
        % (This is the case because the MU is measuring.)
        if strcmp(freq,'24')
            power_from_RU_at_refl = RU(nn).Pt + RU(nn).PL024_LOS - 10*RU(nn).n24_LOS*log10(dist_RU_refl);
            dist_refl_MU = d0*10.^((RSS(1,nn) - power_from_RU_at_refl)./(-10*RU(nn).n24_NLOS));
        elseif strcmp(freq,'58')
            power_from_RU_at_refl = RU(nn).Pt + RU(nn).PL058_LOS - 10*RU(nn).n58_LOS*log10(dist_RU_refl);
            dist_refl_MU = d0*10.^((RSS(1,nn) - power_from_RU_at_refl)./(-10*RU(nn).n58_NLOS));
        end

        % Calculate the total distance:
        total_dist_from_RU2refl2MU(nn) = dist_RU_refl + dist_refl_MU;
        dist_RSS(1,nn) = total_dist_from_RU2refl2MU(nn);
        
    end
end

return
