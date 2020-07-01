function CF = return_connection_Gutman(x,k,l)

phi1=x(1);
phi2=x(2);
theta=x(5);

alpha0=theta;
alpha1=theta+phi1;
alpha2=theta+phi2;

T0=eye(3);
E0=zeros(3,2);
T1=[1,  0,  -0.5*l*sin(alpha0)-0.5*l*sin(alpha1);
    0,  1,  0.5*l*cos(alpha0)+0.5*l*sin(alpha1);
    0,  0,  1,];

E1=[-0.5*l*sin(alpha1),     0;
    0.5*l*cos(alpha1),      0;
    1,                      0];

T1=[1,  0,  -0.5*l*sin(alpha0)-0.5*l*sin(alpha1);
    0,  1,  0.5*l*cos(alpha0)+0.5*l*sin(alpha1);
    0,  0,  1,];

T2=[1,  0,  0.5*l*sin(alpha0)+0.5*l*sin(alpha1);
    0,  1,  -0.5*l*cos(alpha0)-0.5*l*sin(alpha1);
    0,  0,  1,];

E2=[0,      -0.5*l*sin(alpha1);
    0,      0.5*l*cos(alpha1);
    0,      -1];

R0=k*l*[1+(sin(alpha0))^2,           -cos(alpha0)*sin(alpha0),   0;
        -cos(alpha0)*sin(alpha0),   1+(cos(alpha0))^2,          0;
        0,                          0,                          1/6*l^2];
R1=k*l*[1+(sin(alpha1))^2,           -cos(alpha1)*sin(alpha1),   0;
        -cos(alpha1)*sin(alpha1),   1+(cos(alpha1))^2,           0;
        0,                          0,                          1/6*l^2];
R2=k*l*[1+(sin(alpha2))^2,           -cos(alpha2)*sin(alpha2),   0;
        -cos(alpha2)*sin(alpha2),   1+(cos(alpha2))^2,          0;
        0,                          0,                          1/6*l^2];

omega1=T0'*R0*T0+T1'*R1*T1+T2'*R2*T2;
omega2=T0'*R0*E0+T1'*R1*E1+T2'*R2*E2;

CF=-inv(omega1)*omega2;
