function [value, isterminal, direction] = hill4bpSectionEvent(~, s, params)
%HILL4BPSECTIONEVENT  Detect x2=0 crossings with positive x2dot.

r = max(s(1), 0);
theta = s(2);

value = r^params.gamma*sin(theta);
isterminal = 0;
direction = 1;

end
