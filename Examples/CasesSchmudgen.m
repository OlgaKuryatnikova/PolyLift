% these cases are made to run BENCHMARKS for the Schmudgen examples, in SOSTOOLS
% all UPPER BOUND constraints present in the initial cases are INCLUDED in them, and lower bound constraints are non-negativity in those cases; 
% some cases are scaled for numerical stability


%% Example 1
Exper{1}.name = 'ex1';
Exper{1}.numVars=5;
x = varsVector('x',Exper{1}.numVars); 
Exper{1}.f = 7*(2*x(1)-x(2)+x(3)-2*x(4)-2*x(5));
Exper{1}.g = [(7*x(1)-2)^2-49*x(2)^2-(7*x(3)-1)^2-(7*x(5)-1)^2; 49*x(1)*x(3)-49*x(4)*x(5)+49*x(1)^2-1;7*x(3)-49*x(2)^2-49*x(4)^2-1;...
        49*x(1)*x(5)-49*x(2)*x(3)-2;2-sum(x);transpose(x)];
Exper{1}.M = 2;
Exper{1}.L = zeros(Exper{1}.numVars,1);
Exper{1}.U = Exper{1}.M*ones(Exper{1}.numVars,1);
Exper{1}.h = [];

%%  Example 2
Exper{2}.name = 'ex2';
Exper{2}.numVars=10;
x = varsVector('x',Exper{2}.numVars); 
Exper{2}.f = 2*(-x(1)-x(2)+x(3)-2*x(4)-x(5)+x(6)+x(7)-x(8)+x(9)-2*x(10));
Exper{2}.g = [(2*x(3)-2)^2-(2*x(5)-1)^2-2*2*x(6)+4*x(8)^2-(2*x(9)-2)^2+4;-4*x(2)^2+4*x(3)*x(10)-4*x(4)^2+4*x(6)*x(7)-1;...
        4*x(1)*x(8)-4*x(2)*x(3)+4*x(4)*x(7)-4*x(5)*x(10)-2;2.5-sum(x);transpose(x)];
Exper{2}.M = 2.5;
Exper{2}.L = zeros(Exper{2}.numVars,1);
Exper{2}.U = Exper{2}.M*ones(Exper{2}.numVars,1);
Exper{2}.h = [];

%% Example 3
Exper{3}.name = 'ex3';
Exper{3}.numVars=15;
x = varsVector('x',Exper{3}.numVars); 
Exper{3}.f =5*(x(1)-x(2)+x(3)-x(4)-x(5)+x(6)+x(7)-x(8)+x(9)-x(10)+x(11)-x(12)+x(13)-x(14)+x(15)); 
Exper{3}.g = [(5*x(1)-2)^2-25*x(2)^2+(5*x(3)-2)^2-(5*x(4)-1)^2-(5*x(5)-1)^2+(5*x(6)-2)^2-(5*x(7)-2)^2-25*x(8)^2-(5*x(9)-2)^2-(5*x(10)-1)^2+25*x(11)^2-25*x(12)^2+(5*x(13)-2)^2+25*x(14)^2-(5*x(15)-1)^2;...
        -25*x(1)*x(7)-25*x(13)^2-25*x(4)*x(5)+25*x(6)*x(9)+25*x(10)*x(12)-3;25*x(2)*x(3)-25*x(8)*x(11)+25*x(5)*x(15)-25*x(14)^2-3;2-sum(x);transpose(x)];
Exper{3}.M = 2;
Exper{3}.L = zeros(Exper{3}.numVars,1);
Exper{3}.U = Exper{3}.M*ones(Exper{3}.numVars,1);
Exper{3}.h = [];

%% Examples from http://www.minlplib.org/

%% ex_4_1_9
Exper{4}.name = '4_1_9';
Exper{4}.numVars=2;
x = varsVector('x',Exper{4}.numVars); 
Exper{4}.f=-x(1)-x(2);
Exper{4}.g=[-(x(2)-8*x(1)^2+8*x(1)^3-2*x(1)^4)+2;-(x(2)+96*x(1)-88*x(1)^2+32*x(1)^3-4*x(1)^4)+36;3-x(1);4-x(2);transpose(x)];
    %;x(1)-2.3295;x(2)-3.17846];
Exper{4}.M = 7;
Exper{4}.L = zeros(Exper{4}.numVars,1);
Exper{4}.U = [3;4];
Exper{4}.h = [];

%% ex_3_1_4, divided x by their initial bounds
Exper{5}.name = '3_1_4';
Exper{5}.numVars=3;
x = varsVector('x',Exper{5}.numVars); 
Exper{5}.f=(-4*x(1)+4*x(2)-3*x(3)); % use +5 to get a positive value
Exper{5}.g=[(2*4*x(1)*(2*x(1)-4*x(2)+3*x(3))+4*2*x(2)*(4*x(2)-3*x(3))+9*2*x(3)^2+(-2*20*x(1)+4*9*x(2)-3*13*x(3))+24)/4;(-(2*x(1)+4*x(2)+3*x(3))+4)/1;-(4*x(2)+x(3))+2];
Exper{5}.g(4)=1-x(1); % 2
Exper{5}.g(5)=1-x(3); % 3
Exper{5}.g = [Exper{5}.g;transpose(x)];
Exper{5}.M = 3;
Exper{5}.L = zeros(Exper{5}.numVars,1);
Exper{5}.U = [1;1;1];
Exper{5}.h = [];
% the above polynomials are scaled to work for the upper bound hierarchy

%% ex2_1_1
% dvide x by 5 and change all polynomials accordingly
Exper{6}.name = '2_1_1';
Exper{6}.numVars=5;
x = varsVector('x',Exper{6}.numVars); 
Exper{6}.f=(42*x(1)+44*x(2)+45*x(3)+47*x(4)+47.5*x(5))+(-100*x(1)^2-100*x(2)^2-100*x(3)^2-100*x(4)^2-100*x(5)^2)/2;
Exper{6}.g=[(-(20*x(1)+12*x(2)+11*x(3)+7*x(4)+4*x(5))+40)/4;1-x(1);1-x(2)];
Exper{6}.g(4)=1-x(3);
Exper{6}.g(5)=1-x(4);
Exper{6}.g(6)=1-x(5);
Exper{6}.g = [Exper{6}.g;transpose(x)];
Exper{6}.M = Exper{6}.numVars;
Exper{6}.L = zeros(Exper{6}.numVars,1);
Exper{6}.U = 1*ones(Exper{6}.numVars,1);
Exper{6}.h = [];

%% ex2_1_2, scaled x(6) -> x(6)/10
% optimum at degree 4, small sos are not redundant
Exper{7}.name = '2_1_2';
Exper{7}.numVars=6;
x = varsVector('x',Exper{7}.numVars); 
Exper{7}.f=-10.5*x(1)-7.5*x(2)-3.5*x(3)-2.5*x(4)-1.5*x(5)-100*x(6)+(-x(1)^2-x(2)^2-x(3)^2-x(4)^2-x(5)^2)/2;
Exper{7}.g=[-(6*x(1)+3*x(2)+3*x(3)+2*x(4)+x(5))+6.5;-(x(1)+x(3)+x(6))+2];
Exper{7}.g(3)=1-x(1);
Exper{7}.g(4)=1-x(2);
Exper{7}.g(5)=1-x(3);
Exper{7}.g(6)=1-x(4);
Exper{7}.g(7)=1-x(5);
Exper{7}.g = [Exper{7}.g;transpose(x)];
Exper{7}.M = 7;
Exper{7}.L = zeros(Exper{7}.numVars,1);
Exper{7}.U = 1*ones(Exper{7}.numVars,1);
Exper{7}.U(end) = Exper{7}.U(end)+1;
Exper{7}.h = [];

%% ex2_1_4
Exper{8}.name = '2_1_4';
Exper{8}.numVars=6;
x = varsVector('x',Exper{8}.numVars); 
Exper{8}.f=6.5*x(1)-x(2)-2*x(3)-3*x(4)-2*x(5)-x(6)-0.5*x(1)^2; 
Exper{8}.g=[-(x(1)+2*x(2)+8*x(3)+x(4)+3*x(5)+5*x(6))+16;-(-8*x(1)-4*x(2)-2*x(3)+2*x(4)+4*x(5)-x(6))-1;-(2*x(1)+0.5*x(2)+0.2*x(3)-3*x(4)-x(5)-4*x(6))+24];
Exper{8}.g(4)=-(0.2*x(1)+2*x(2)+0.1*x(3)-4*x(4)+2*x(5)+2*x(6))+12;
Exper{8}.g(5)=-(-0.1*x(1)-0.5*x(2)+2*x(3)+5*x(4)-5*x(5)+3*x(6))+3;
Exper{8}.g(6)=1-x(1);
Exper{8}.g(7)=1-x(4);
Exper{8}.g(8)=1-x(5);
Exper{8}.g(9)=2-x(6);
Exper{8}.g = [Exper{8}.g;transpose(x)];
Exper{8}.M = 15;
Exper{8}.L = zeros(Exper{8}.numVars,1);
Exper{8}.U = [1;8;2;1;1;2];
Exper{8}.h = [];

%% ex2_1_3
% optimum at degree 4, small sos are not redundant
Exper{9}.name = '2_1_3';
Exper{9}.numVars=13;
x = varsVector('x',Exper{9}.numVars); 
Exper{9}.f=5*x(1)+5*x(2)+5*x(3)+5*x(4)-x(5)-x(6)-x(7)-x(8)-x(9)-x(10)-x(11)-x(12)-x(13)+...
    (-10*x(1)^2-10*x(2)^2-10*x(3)^2-10*x(4)^2)/2;
Exper{9}.g = [-(2*x(1)+2*x(2)+x(10)+x(11))+10;-(2*x(1)+2*x(3)+x(10)+x(12))+10;-(2*x(2)+2*x(3)+x(11)+x(12))+10];
Exper{9}.g(4)=-(-8*x(1)+x(10));
Exper{9}.g(5)=-(-8*x(2)+x(11));
Exper{9}.g(6)=-(-8*x(3)+x(12));
Exper{9}.g(7)=-(-2*x(4)-x(5)+x(10));
Exper{9}.g(8)=-(-2*x(6)-x(7)+x(11));
Exper{9}.g(9)=-(-2*x(8)-x(9)+x(12));
Exper{9}.g(10)=1-x(1);
Exper{9}.g(11)=1-x(2);
Exper{9}.g(12)=1-x(3);
Exper{9}.g(13)=1-x(4);
Exper{9}.g(14)=1-x(5);
Exper{9}.g(15)=1-x(6);
Exper{9}.g(16)=1-x(7);
Exper{9}.g(17)=1-x(8);
Exper{9}.g(18)=1-x(9);
Exper{9}.g(19)=1-x(13);
Exper{9}.g = [Exper{9}.g;transpose(x)];
Exper{9}.M = 19;
Exper{9}.L = zeros(Exper{9}.numVars,1);
Exper{9}.U = 1*ones(Exper{9}.numVars,1);
Exper{9}.U(end-2:end) = Exper{9}.U(end-2:end)+2;
Exper{9}.h = [];

%% ex2_1_7
Exper{10}.name = '2_1_7';
Exper{10}.numVars=20;
x = varsVector('x',Exper{10}.numVars); 
Exper{10}.f=2*x(1)+4*x(2)+6*x(3)+8*x(4)+10*x(5)+12*x(6)+14*x(7)+16*x(8)+18*x(9)+20*x(10)+22*x(11)+24*x(12)+26*x(13)+...
    28*x(14)+30*x(15)+32*x(16)+34*x(17)+36*x(18)+38*x(19)+40*x(20)+(-x(1)^2-2*x(2)^2-3*x(3)^2-4*x(4)^2-...
    5*x(5)^2-6*x(6)^2-7*x(7)^2-8*x(8)^2-9*x(9)^2-10*x(10)^2-11*x(11)^2-12*x(12)^2-13*x(13)^2-14*x(14)^2-...
    15*x(15)^2-16*x(16)^2-17*x(17)^2-18*x(18)^2-19*x(19)^2-20*x(20)^2)/2-4.2e2;
Exper{10}.g=[-(-3*x(1)+7*x(2)-5*x(4)+x(5)+x(6)+2*x(8)-x(9)-x(10)-9*x(11)+3*x(12)+5*x(13)+x(16)+7*x(17)-7*x(18)-4*x(19)-6*x(20))-5;...
    -(7*x(1)-5*x(3)+x(4)+x(5)+2*x(7)-x(8)-x(9)-9*x(10)+3*x(11)+5*x(12)+x(15)+7*x(16)-7*x(17)-4*x(18)-6*x(19)-3*x(20))+2];
Exper{10}.g(3)=-(-5*x(2)+x(3)+x(4)+2*x(6)-x(7)-x(8)-9*x(9)+3*x(10)+5*x(11)+x(14)+7*x(15)-7*x(16)-4*x(17)-6*x(18)-3*x(19)+7*x(20))-1;
Exper{10}.g(4)=-(-5*x(1)+x(2)+x(3)+2*x(5)-x(6)-x(7)-9*x(8)+3*x(9)+5*x(10)+x(13)+7*x(14)-7*x(15)-4*x(16)-6*x(17)-3*x(18)+7*x(19))-3;
Exper{10}.g(5)=-(x(1)+x(2)+2*x(4)-x(5)-x(6)-9*x(7)+3*x(8)+5*x(9)+x(12)+7*x(13)-7*x(14)-4*x(15)-6*x(16)-3*x(17)+7*x(18)-5*x(20))+5;
Exper{10}.g(6)=-(x(1)+2*x(3)-x(4)-x(5)-9*x(6)+3*x(7)+5*x(8)+x(11)+7*x(12)-7*x(13)-4*x(14)-6*x(15)-3*x(16)+7*x(17)-5*x(19)+x(20))+4;
Exper{10}.g(7)=-(2*x(2)-x(3)-x(4)-9*x(5)+3*x(6)+5*x(7)+x(10)+7*x(11)-7*x(12)-4*x(13)-6*x(14)-3*x(15)+7*x(16)-5*x(18)+x(19)+x(20))-1;
Exper{10}.g(8)=-(2*x(1)-x(2)-x(3)-9*x(4)+3*x(5)+5*x(6)+x(9)+7*x(10)-7*x(11)-4*x(12)-6*x(13)-3*x(14)+7*x(15)-5*x(17)+x(18)+x(19));
Exper{10}.g(9)=-(-x(1)-x(2)-9*x(3)+3*x(4)+5*x(5)+x(8)+7*x(9)-7*x(10)-4*x(11)-6*x(12)-3*x(13)+7*x(14)-5*x(16)+x(17)+x(18)+2*x(20))+9;
Exper{10}.g(10)=-(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)+x(11)+x(12)+x(13)+x(14)+x(15)+x(16)+x(17)+x(18)+x(19)+x(20))+40;
Exper{10}.g = [Exper{10}.g;transpose(x)];
Exper{10}.h = [];
Exper{10}.M = 40;
Exper{10}.L = zeros(Exper{10}.numVars,1);
Exper{10}.U = 40*ones(Exper{10}.numVars,1);

%% ex2_1_5
Exper{11}.name = '2_1_5';
Exper{11}.numVars=10;
x = varsVector('x',Exper{11}.numVars); 
Exper{11}.f = -20*x(1)-80*x(2)-20*x(3)-50*x(4)-60*x(5)-90*x(6)+10*x(8)+10*x(9)+10*x(10)-5*x(1)^2-5*x(2)^2-5*x(3)^2-5*x(4)^2-5*x(5)^2-5*x(6)^2-5*x(7)^2;
Exper{11}.g=[-(-2*x(1)-6*x(2)-x(3)-3*x(5)-3*x(6)-2*x(7)-6*x(8)-2*x(9)-2*x(10))-4;-(6*x(1)-5*x(2)+8*x(3)-3*x(4)+x(6)+3*x(7)+8*x(8)+9*x(9)-3*x(10))+22];
Exper{11}.g(3)=-(-5*x(1)+6*x(2)+5*x(3)+3*x(4)+8*x(5)-8*x(6)+9*x(7)+2*x(8)-9*x(10))-6;
Exper{11}.g(4)=-(9*x(1)+5*x(2)-9*x(4)+x(5)-8*x(6)+3*x(7)-9*x(8)-9*x(9)-3*x(10))-23;
Exper{11}.g(5)=-(-8*x(1)+7*x(2)-4*x(3)-5*x(4)-9*x(5)+x(6)-7*x(7)-x(8)+3*x(9)-2*x(10))-12;
Exper{11}.g(6)=-(-7*x(1)-5*x(2)-2*x(3)-6*x(5)-6*x(6)-7*x(7)-6*x(8)+7*x(9)+7*x(10))-3;
Exper{11}.g(7)=-(x(1)-3*x(2)-3*x(3)-4*x(4)-x(5)-4*x(7)+x(8)+6*x(9))+1;
Exper{11}.g(8)=-(x(1)-2*x(2)+6*x(3)+9*x(4)-7*x(6)+9*x(7)-9*x(8)-6*x(9)+4*x(10))+12;
Exper{11}.g(9)=-(-4*x(1)+6*x(2)+7*x(3)+2*x(4)+2*x(5)+6*x(7)+6*x(8)-7*x(9)+4*x(10))+15;
Exper{11}.g(10)=-(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10))+9;
Exper{11}.g(11)=-(-x(1)-x(2)-x(3)-x(4)-x(5)-x(6)-x(7)-x(8)-x(9)-x(10))-1;
Exper{11}.g(12)=1-x(1);
Exper{11}.g(13)=1-x(2);
Exper{11}.g(14)=1-x(3);
Exper{11}.g(15)=1-x(4);
Exper{11}.g(16)=1-x(5);
Exper{11}.g(17)=1-x(6);
Exper{11}.g(18)=1-x(7);
Exper{11}.g(19)=1-x(8);
Exper{11}.g(20)=1-x(9);
Exper{11}.g(21)=1-x(10);
Exper{11}.g = [Exper{11}.g;transpose(x)];
Exper{11}.h = [];
Exper{11}.M = 1;
Exper{11}.L = zeros(Exper{11}.numVars,1);
Exper{11}.U = 1*ones(Exper{11}.numVars,1);

%% ex2_1_6
Exper{12}.name = '2_1_6';
Exper{12}.numVars=10;
x = varsVector('x',Exper{12}.numVars); 
Exper{12}.f =48*x(1)+42*x(2)+48*x(3)+45*x(4)+44*x(5)+41*x(6)+47*x(7)+42*x(8)+45*x(9)+46*x(10)-50*x(1)^2-50*x(2)^2-50*x(3)^2- ...
    50*x(4)^2-50*x(5)^2-50*x(6)^2-50*x(7)^2-50*x(8)^2-50*x(9)^2-50*x(10)^2;
Exper{12}.g=[-(- 2*x(1)-6*x(2)-x(3)-3*x(5)-3*x(6)-2*x(7)-6*x(8)-2*x(9)-2*x(10))-4;...
-(6*x(1)-5*x(2)+8*x(3)-3*x(4)+x(6)+3*x(7)+8*x(8)+9*x(9)-3*x(10))+22;...
-(-5*x(1)+6*x(2)+5*x(3)+3*x(4)+8*x(5)-8*x(6)+9*x(7)+2*x(8)-9*x(10))-6;...
-(9*x(1)+5*x(2)-9*x(4)+x(5)-8*x(6)+3*x(7)-9*x(8)-9*x(9)-3*x(10))-23;...
-(-8*x(1)+7*x(2)-4*x(3)-5*x(4)-9*x(5)+x(6)-7*x(7)-x(8)+3*x(9)-2*x(10))-12];
Exper{12}.g(6)=1-x(1);
Exper{12}.g(7)=1-x(2);
Exper{12}.g(8)=1-x(3);
Exper{12}.g(9)=1-x(4);
Exper{12}.g(10)=1-x(5);
Exper{12}.g(11)=1-x(6);
Exper{12}.g(12)=1-x(7);
Exper{12}.g(13)=1-x(8);
Exper{12}.g(14)=1-x(9);
Exper{12}.g(15)=1-x(10);
Exper{12}.g = [Exper{12}.g;transpose(x)];
Exper{12}.h=[];
Exper{12}.M = 10;
Exper{12}.L = zeros(Exper{12}.numVars,1);
Exper{12}.U = 1*ones(Exper{12}.numVars,1);

%% ex3_1_3
Exper{13}.name = '3_1_3';
Exper{13}.numVars=6;
x = varsVector('x',Exper{13}.numVars); 
Exper{13}.f =-138+100*x(1)+4*x(2)+2*x(3)+8*x(4)+2*x(5)+8*x(6)-25*x(1)^2-x(2)^2-x(3)^2-x(4)^2-x(5)^2-x(6)^2;
Exper{13}.g =[-6*x(3)+x(4)+x(3)^2+5;-6*x(5)+x(6)+x(5)^2+5;-(x(1)-3*x(2))+2;-(-x(1)+x(2))+2;-(x(1)+x(2))+6;x(1)+x(2)-2];
Exper{13}.g = [Exper{13}.g;x(1);x(2);x(3)-1;5-x(3);6-x(4);x(4);x(5)-1;5-x(5);10-x(6);x(6)];
Exper{13}.L = [0;0;1;0;1;0];
Exper{13}.U = [6;4;5;6;5;10];
Exper{13}.M = sum(Exper{13}.U);
Exper{13}.h=[];

%% mathopt2
Exper{14}.name = 'mathopt2';
Exper{14}.numVars=2;
x = varsVector('x',Exper{14}.numVars); 
Exper{14}.f = (40*(4*x(1)-3)^4-16*(4*x(1)-3)^2*(2*x(2)-1)+2*(2*x(2)-1)^2)/16;
Exper{14}.h = [(4*x(1)-3)-10*(2*x(2)-1)-(4*x(1)-3)*(2*x(2)-1);(4*x(1)-3)-3*(2*x(2)-1)];
Exper{14}.g = [(-((4*x(1)-3)+(2*x(2)-1))+1)/4; ((4*x(1)-3)-(2*x(2)-1)+2)/4];
Exper{14}.M = 2;
Exper{14}.L = zeros(Exper{14}.numVars,1);
Exper{14}.U = [1;1]; % upper bounds are lose but work well numerically
Exper{14}.M = 2;
Exper{14}.L = zeros(Exper{14}.numVars,1);
Exper{14}.U = [1;1];

%% ex2_1_10
Exper{15}.name = 'ex2_1_10';
Exper{15}.numVars=20;
x = varsVector('x',Exper{15}.numVars); 
Exper{15}.f=547663.5+100*(-1197*x(1)-405*x(2)-1012*x(3)-4823*x(4)-1.89e3*x(5)+1.3e3*x(6)-2937*x(7)-1334*x(8)+3526*x(9)+...
    1558*x(10)+2184*x(11)+294*x(12)-3888*x(13)-2.73e3*x(14)+935*x(15)-4284*x(16)-1647*x(17)+4941*x(18)-...
    3686*x(19)+1898*x(20))+10000*(-31.5*x(1)^2-7.5*x(2)^2-22*x(3)^2-45.5*x(4)^2-22.5*x(5)^2-25*x(6)^2-44.5*x(7)^2-...
29*x(8)^2-43*x(9)^2-41*x(10)^2+21*x(11)^2+49*x(12)^2+24*x(13)^2+45.5*x(14)^2+5.5*x(15)^2+31.5*x(16)^2+...
30.5*x(17)^2+30.5*x(18)^2+19*x(19)^2+13*x(20)^2); 
Exper{15}.f = Exper{15}.f/10000;
Exper{15}.g=[-100*(3*x(1)+5*x(2)+5*x(3)+6*x(4)+4*x(5)+4*x(6)+5*x(7)+6*x(8)+4*x(9)+4*x(10)+8*x(11)+4*x(12)+2*x(13)+x(14)+x(15)+x(16)+2*x(17)+x(18)+7*x(19)+3*x(20))+3.8e2;...
 -100*(5*x(1)+4*x(2)+5*x(3)+4*x(4)+x(5)+4*x(6)+4*x(7)+2*x(8)+5*x(9)+2*x(10)+3*x(11)+6*x(12)+x(13)+7*x(14)+7*x(15)+5*x(16)+8*x(17)+7*x(18)+2*x(19)+x(20))+415;...
 -100*(x(1)+5*x(2)+2*x(3)+4*x(4)+7*x(5)+3*x(6)+x(7)+5*x(8)+7*x(9)+6*x(10)+x(11)+7*x(12)+2*x(13)+4*x(14)+7*x(15)+5*x(16)+3*x(17)+4*x(18)+x(19)+2*x(20))+385;...
 -100*(3*x(1)+2*x(2)+6*x(3)+3*x(4)+2*x(5)+x(6)+6*x(7)+x(8)+7*x(9)+3*x(10)+7*x(11)+7*x(12)+8*x(13)+2*x(14)+3*x(15)+4*x(16)+5*x(17)+8*x(18)+x(19)+2*x(20))+405;...
 -100*(6*x(1)+6*x(2)+6*x(3)+4*x(4)+5*x(5)+2*x(6)+2*x(7)+4*x(8)+3*x(9)+2*x(10)+7*x(11)+5*x(12)+3*x(13)+6*x(14)+7*x(15)+5*x(16)+8*x(17)+4*x(18)+6*x(19)+3*x(20))+4.7e2;...
 -100*(5*x(1)+5*x(2)+2*x(3)+x(4)+3*x(5)+5*x(6)+5*x(7)+7*x(8)+4*x(9)+3*x(10)+4*x(11)+x(12)+7*x(13)+3*x(14)+8*x(15)+3*x(16)+x(17)+6*x(18)+2*x(19)+8*x(20))+415;...
 -100*(3*x(1)+6*x(2)+6*x(3)+3*x(4)+x(5)+6*x(6)+x(7)+6*x(8)+7*x(9)+x(10)+4*x(11)+3*x(12)+x(13)+4*x(14)+3*x(15)+6*x(16)+4*x(17)+6*x(18)+5*x(19)+4*x(20))+400;...
 -100*(x(1)+2*x(2)+x(3)+7*x(4)+8*x(5)+7*x(6)+6*x(7)+5*x(8)+8*x(9)+7*x(10)+2*x(11)+3*x(12)+5*x(13)+5*x(14)+4*x(15)+5*x(16)+4*x(17)+2*x(18)+2*x(19)+8*x(20))+4.6e2;...
 -100*(8*x(1)+5*x(2)+2*x(3)+5*x(4)+3*x(5)+8*x(6)+x(7)+3*x(8)+3*x(9)+5*x(10)+4*x(11)+5*x(12)+5*x(13)+6*x(14)+x(15)+7*x(16)+x(17)+2*x(18)+2*x(19)+4*x(20))+400];
Exper{15}.g=[Exper{15}.g/100;-(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)+x(11)+x(12)+x(13)+x(14)+x(15)+x(16)+x(17)+x(18)+x(19)+x(20))+2;transpose(x)];
Exper{15}.M = 2;
Exper{15}.L = zeros(Exper{15}.numVars,1);
Exper{15}.U = 2*ones(Exper{15}.numVars,1);
Exper{15}.h=[];