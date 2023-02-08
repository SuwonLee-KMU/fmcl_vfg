classdef VFGParametric
    % Vector field based guidance algorithm for parametric curve in space.

    properties
        curveObj
        conv_factor
    end

    methods
        % Constructor method
        function obj = VFGParametric(curveObj, conv_factor)
            obj.curveObj    = curveObj;
            obj.conv_factor = conv_factor;
        end

        % Evaluation method
        function v_d = feval(obj, pos_vector)
            sz = size(pos_vector);
            if sz(2) ~= 3
                error('pos_vector should be in shape of [n x 3]!')
            end
            npoints = sz(1);
            v_d     = zeros(npoints,3);
            for i = 1:npoints
                pos = pos_vector(i,:);
                nvc = obj.feval_for_one_vector(pos);
                v_d(i,:) = nvc;
            end
        end

        % Evaluation for one point
        function v_d = feval_for_one_vector(obj, pos_vector)
            [pos_near, tau] = obj.nearest_param_from(pos_vector);
            basis           = obj.basis_vectors_at(tau, pos_vector);
            coef            = obj.dist_based_coef(norm(pos_near-pos_vector));
            vel_vector      = coef(1)*basis(1,:) + coef(2)*basis(2,:);
            v_d             = vel_vector/norm(vel_vector);
        end
        
        % Compute cofficient based on distance to the curve
        function coef = dist_based_coef(obj, dist)
            k_conv = atan(dist*obj.conv_factor)/pi*2;
            k_trav = sqrt(1-k_conv^2);
            coef   = [k_conv, k_trav];
        end
        
        % Calculate nearest point on the curve from a space position.
        function [nearest_point, tau_value] = nearest_param_from(obj, pos)
            fun               = @(tau) norm(obj.curveObj.feval(tau) - pos);
            tau_range         = [0,1];
            [tau_value, fval] = fminbnd(fun, tau_range(1), tau_range(2));
            nearest_point     = obj.curveObj.feval(tau_value);
        end

        % Compute basis vector at curve point of param=tau with respect to a space position.
        function basis = basis_vectors_at(obj, tau, pos)
            if numel(tau) > 1 
                error('tau should be a scalar.');
            end
            tan = obj.curveObj.tangent(tau);
            tan = tan/norm(tan);
            cpos = obj.curveObj.feval(tau);
            basis1 = (cpos - pos)/norm(cpos-pos);
            basis2 = tan;
            basis = [basis1; basis2];
        end

        % for visualization
        function positions = generate_random_position_near_curve(obj, npoints, scale)
            ul = 1;
            ll = 0;
            rand_taus = rand(npoints, 1)*(ul-ll)+ll;
            positions = zeros(npoints, 3);
            for i = 1:npoints
                tau            = rand_taus(i);
                pos_on_curve   = obj.curveObj.feval(tau);
                random_pos     = pos_on_curve + randn(1,3)*scale;
                positions(i,:) = random_pos;
            end
        end
    end
end