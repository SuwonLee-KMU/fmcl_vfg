classdef SpaceCurvePoly
    % Class for space curve.

    properties
        coef_x
        coef_y
        coef_z
    end

    methods
        % Constructor method
        function obj = SpaceCurvePoly(coef_x, coef_y, coef_z)
            obj.coef_x = coef_x;
            obj.coef_y = coef_y;
            obj.coef_z = coef_z;
        end
       
        % Evaluation method
        function pos = feval(obj, tau)
            if numel(tau) == 1
                x = polyval(obj.coef_x,tau);
                y = polyval(obj.coef_y,tau);
                z = polyval(obj.coef_z,tau);
                pos = [x, y, z];
            else
                pos = zeros(numel(tau),3);
                for i = 1:numel(tau)
                    t = tau(i);
                    x = polyval(obj.coef_x,t);
                    y = polyval(obj.coef_y,t);
                    z = polyval(obj.coef_z,t);
                    pos(i,:) = [x, y, z];
                end
            end
        end

        % Normalized tangent vector at given point
        function tan = tangent(obj, tau)
            if numel(tau) == 1
                tx = polyval(polyder(obj.coef_x), tau);
                ty = polyval(polyder(obj.coef_y), tau);
                tz = polyval(polyder(obj.coef_z), tau);
                tan = [tx, ty, tz];
            else
                tan = zeros(numel(tau),3);
                for i = 1:numel(tau)
                    t = tau(i);
                    x = polyval(polyder(obj.coef_x),t);
                    y = polyval(polyder(obj.coef_y),t);
                    z = polyval(polyder(obj.coef_z),t);
                    tan(i,:) = [x, y, z];
                end
            end
        end
    end
end


    