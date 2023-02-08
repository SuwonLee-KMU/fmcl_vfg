% Author: Suwon Lee from Kookmin University
% Date: February 7, 2023 6:01:06 PM GMT+9
% Code: Class for multiple curve objects.

classdef MultipleCurvesPoly
    properties
        curves
    end

    methods
        function obj = MultipleCurvesPoly(curves_cell) % Constructor method
            obj.curves = curves_cell;
        end

        function pos = feval(obj, tau) % Evaluation method
            [tau, idx] = obj.rescale_param(tau);
            pos = zeros(numel(tau),3);
            for i = 1:numel(tau)
                pos(i,:) = obj.curves{idx(i)}.feval(tau(i));
            end
        end

        function pos = feval_for_one_tau(obj, tau) % Evaluation method for one point
            [tau, idx] = obj.rescale_param(tau);
            pos = obj.curves{idx}.feval(tau);
        end

        function tan = tangent(obj, tau) % Normalized tangent vector at given point.
            [tau, idx] = obj.rescale_param(tau);
            tan = zeros(numel(tau),3);
            for i = 1:numel(tau)
                tan(i,:) = obj.curves{idx(i)}.tangent(tau(i));
            end
        end

        function [tau_new, curve_index] = rescale_param(obj, tau) % Rescale parameter for multiple curve into [0, 1]
            tau_new    = zeros(numel(tau),1);
            cuve_index = zeros(numel(tau),1);
            for i = 1:numel(tau)
                [tn, ci]       = obj.rescale_param_for_one_tau(tau(i));
                tau_new(i)     = tn;
                curve_index(i) = ci;
            end
        end

        function [tau_new, curve_index] = rescale_param_for_one_tau(obj, tau) % Rescale parameter for multiple curve into [0, 1]
            ncurves     = numel(obj.curves);
            expanded    = tau*ncurves;
            curve_index = max(1,ceil(expanded));
            tau_new     = expanded - (curve_index - 1);
        end
    end
end