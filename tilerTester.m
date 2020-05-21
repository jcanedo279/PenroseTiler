n = 5;
k = 10;
shiftVect = [1E-10 1E-10 1E-10 1E-10 1E-10];
%shiftVect = normShiftVect(0, sv);

setGlobals(n, k, shiftVect);

genGrid();


function setGlobals(nVal, kVal, shiftVectVal)
    global n
    global k
    global shiftVect

    % Specify Multigrid Dimmensions
    n = nVal;
    k = kVal;

    % Create shift vector of shifts in each of the five directions
    shiftVect = shiftVectVal;
end

function genGrid()
    global n
	global k
    
    figure(1)
    for r = 0 : n-1
        for a = -k : k
            for s = r+1 : n-1
                for b = -k : k
                    % p_rs_ab
                    p_rs_ab_vect = pointCord(r, s, a, b);
                    p_rs_ab_tileVerts = tileCord(r, s, p_rs_ab_vect);
                    p_rs_ab_spaceCordVect = spaceCordVect(p_rs_ab_vect);
                    if (p_rs_ab_spaceCordVect(r+1) ~= a) || (p_rs_ab_spaceCordVect(s+1) ~= b)
                        disp([r s a b])
                        disp(p_rs_ab_spaceCordVect())
                    end
                    
                    if b ~= 0
                        % p_rs_a-b
                        p_rs_anb_vect = pointCord(r, s, a, -b);
                        p_rs_anb_tileVerts = tileCord(r, s, p_rs_anb_vect);
                        p_rs_anb_spaceCordVect = spaceCordVect(p_rs_anb_vect);
                        if (p_rs_anb_spaceCordVect(r+1) ~= a) || (p_rs_anb_spaceCordVect(s+1) ~= -b)
                            disp([r s a b])
                            disp(p_rs_ab_spaceCordVect())
                        end
                    
                    elseif a ~= 0
                        % p_rs_-ab
                        p_rs_nab_vect = pointCord(r, s, -a, b);
                        p_rs_nab_tileVerts = tileCord(r, s, p_rs_nab_vect);
                        p_rs_nab_spaceCordVect = spaceCordVect(p_rs_nab_vect);
                        if (p_rs_nab_spaceCordVect(r+1) ~= -a) || (p_rs_nab_spaceCordVect(s+1) ~= b)
                            disp([r s -a b])
                            disp(p_rs_nab_spaceCordVect())
                        end
                    
                    elseif (a ~= 0) && (b ~= 0)
                        % p_rs_-a-b
                        p_rs_nanb_vect = pointCord(r, s, -a, -b);
                        p_rs_nanb_tileVerts = tileCord(r, s, p_rs_nanb_vect);
                        p_rs_nanb_spaceCordVect = spaceCordVect(p_rs_nanb_vect);
                        if (p_rs_nanb_spaceCordVect(r+1) ~= -a) || (p_rs_nanb_spaceCordVect(s+1) ~= -b)
                            disp([r s -a -b])
                            disp(p_rs_nanb_spaceCordVect())
                        end
                    end
                    
                    if a == b
                        patch(p_rs_ab_tileVerts(:,1), p_rs_ab_tileVerts(:,2), 'black');
                    else
                        patch(p_rs_ab_tileVerts(:,1), p_rs_ab_tileVerts(:,2), 'white');
                    end
                        
                    axis equal
                    grid on
                    
                    %keyboard

                end
            end
        end
    end
    
    drawnow
    
end

function eKp = tileCord(r, s, p)
global n

    eKp_x = 0;
    eKp_y = 0;
    for i = 0 : n-1
        nV = normVect(i);
        Ki = spaceCord(p, i);
        eKp_x = eKp_x + nV(1)*Ki;
        eKp_y = eKp_y + nV(2)*Ki;
    end
    eKpc = [eKp_x eKp_y];
    
    nV_r = normVect(r);
    eKp_r_x = eKp_x + nV_r(1);
    eKp_r_y = eKp_y + nV_r(2);
    eKp_r = [eKp_r_x eKp_r_y];
    
    nV_s = normVect(s);
    eKp_s_x = eKp_x + nV_s(1);
    eKp_s_y = eKp_y + nV_s(2);
    eKp_s = [eKp_s_x eKp_s_y];
    
    eKp_rs_x = eKp_x + nV_r(1) + nV_s(1);
    eKp_rs_y = eKp_y + nV_r(2) + nV_s(2);
    eKp_rs = [eKp_rs_x eKp_rs_y];
    
    eKp = [eKpc; eKp_r; eKp_rs; eKp_s; eKpc];
end

function Kp = spaceCordVect(p)
    global n
    Kp = [0 0 0 0 0];
    
    for i = 0 : n-1
        Kp(i+1) = spaceCord(p, i);
    end
end

function Ki = spaceCord(p, i)
    global shiftVect
    
    nV = normVect(i);
    Ki = round(p(1)*nV(1) + p(2)*nV(2) - shiftVect(i+1));
    
end

function p = pointCord(r, s, a, b)
    zr = zeta(r);
    zs = zeta(s);
    
    x_rs_num = phiOp(r, a)*(sin(zs)) - phiOp(s, b)*(sin(zr));
    x_rs_den = cos(zr)*sin(zs) - cos(zs)*sin(zr);
    x_rs = x_rs_num / x_rs_den;
    
    y_rs_num = phiOp(r, a)/cos(zr) - phiOp(s, b)/cos(zs);
    y_rs_den = tan(zr) - tan(zs);
    y_rs = y_rs_num / y_rs_den;
    
    p = [x_rs y_rs];
end

function phi = phiOp(i, t)
    global shiftVect
    phi = t + (1/2) - shiftVect(i+1);
end

function e = normVect(i)
    e_x = cos(zeta(i));
    e_y = sin(zeta(i));
    e = [e_x e_y];
end

function z = zeta(i)
    z = (2/5)*(pi)*i;
end

function sV = normShiftVect(goalSum, vector)
    global n
    currentSum = 0;
    for i = 0 : n-1
        currentSum = currentSum + vector(i+1);
    end
    
    sV = [0 0 0 0 0];
    
    if (currentSum == goalSum)
        scalingFactor = goalSum/currentSum;
 
        for i = 0 : n-1
            sV(i+1) = vector(i+1)*scalingFactor;
        end
    end
end




