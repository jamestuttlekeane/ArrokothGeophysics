clear;

% This MATLAB script calculates the gravity field, and a variety of other
% data products, for (486958) Arrokoth. This script was created by Dr.
% James Tuttle Keane (NASA Jet Propulsion Laboratory, California Institute 
% of Technology, james.t.keane@jpl.nasa.gov)

%-- constants
    d2r = 2*pi/360;
    r2d = 360/2/pi;
    G = 6.674e-11; % gravitational constant, m^3 kg^-2 s^-2
        
    disp(' ');
    
%% IMPORTING GLOBAL MESH
disp(' ');disp('importing global mesh');

%-- importing mesh
    global_mesh = 'MU69-Merged-EPSC';
    
    f=stlread([global_mesh,'.stl']);
    
    disp(['importing file: ',global_mesh,'.stl'])

%-- scaling up the mesh (to meters)
    f.vertices=f.vertices*1000;
    
%-- shrink it down a bit (this is useful if you want to run the code fast)
    %mesh_resolution=0.01; % good for external field
    mesh_resolution=1.0;
    if mesh_resolution~=1
    [faces,vertices]=reducepatch(f.faces,f.vertices,mesh_resolution);
    f.faces=faces;
    f.vertices=vertices;
    end
    disp(['number of faces in mesh: ',num2str(size(f.faces,1))])

%-- performing other mesh calculations for future use
    TR = triangulation(f.faces,f.vertices);
    centers = incenter(TR);
    normals = faceNormal(TR);
    area = areaIsosurface(f.faces,f.vertices);

    n=size(centers,1);
    vert1=f.vertices(f.faces(:,1),:);
    vert2=f.vertices(f.faces(:,2),:);
    vert3=f.vertices(f.faces(:,3),:);

    [faces_a, vertices_a] = freeBoundary(TR);
    
% PRINTING GEOMETRIC PROPERTIES

    xxx=[vert1(:,1); vert2(:,1); vert3(:,1); centers(:,1)];
    yyy=[vert1(:,2); vert2(:,2); vert3(:,2); centers(:,2)];
    zzz=[vert1(:,3); vert2(:,3); vert3(:,3); centers(:,3)];
    
    disp(['x: max = ',num2str(max(xxx)/1000,'%3.2f'),' km, min = ',num2str(min(xxx)/1000,'%3.2f'),' km, delta = ',num2str((max(xxx)-min(xxx))/1000,'%3.2f'),' km']);
    disp(['y: max = ',num2str(max(yyy)/1000,'%3.2f'),' km, min = ',num2str(min(yyy)/1000,'%3.2f'),' km, delta = ',num2str((max(yyy)-min(yyy))/1000,'%3.2f'),' km']);
    disp(['z: max = ',num2str(max(zzz)/1000,'%3.2f'),' km, min = ',num2str(min(zzz)/1000,'%3.2f'),' km, delta = ',num2str((max(zzz)-min(zzz))/1000,'%3.2f'),' km']);

%% IMPORTING OBJECT FOR SCALE;

% These shape files were downloaded from Thingiverse. They are only used to
% provide a fun sense of scale. They can be swapped in with other objects
% for scale, or deleted.

% %-- golden gate bridge
    f_SCALE=stlread('golden_gate_bridge.stl');
    SCALE_LENGTH = max(f_SCALE.vertices(:,1))-min(f_SCALE.vertices(:,1));
    f_SCALE.vertices(:,1)=f_SCALE.vertices(:,1) * 1970/SCALE_LENGTH;

    SCALE_LENGTH = max(f_SCALE.vertices(:,2))-min(f_SCALE.vertices(:,2));
    f_SCALE.vertices(:,2)=f_SCALE.vertices(:,2) * 27/SCALE_LENGTH;

    SCALE_LENGTH = max(f_SCALE.vertices(:,3))-min(f_SCALE.vertices(:,3));
    f_SCALE.vertices(:,3)=f_SCALE.vertices(:,3) * 230/SCALE_LENGTH;
    
%-- enterprise d
    f_SCALE2=stlread('Enterprise-D.stl');
    SCALE_LENGTH = max(f_SCALE2.vertices(:,1))-min(f_SCALE2.vertices(:,1));
    f_SCALE2.vertices=f_SCALE2.vertices * 641/SCALE_LENGTH;

%% QUICK VIEW
    
figure(1);clf;hold on;
set(gcf,'color','w');
fonts=16;
fontn='helvetica neue';
set(gca,'fontsize',fonts,'fontname',fontn);

%-- figure size
    set(gcf,'Units','centimeters')    
    set(gcf,'position',[15 24 40 22.5]);
    set(gcf,'PaperUnits','centimeters');
    P0=get(gcf,'position');
    set(gcf,'PaperPosition',[0 0 P0(3) P0(4);])
    
    %-- lighting properties (stark)
        light_ambient=0.0;           % fractional strength of ambient light
        light_diffuse=0.99;           % fractional strength of diffuse light
        light_specular=0.01;          % fractional strength of specular light
        light_specularexponent=0.1;    % specular light exponent
        
        
    %-- trisurf (shape model)
        t1=trisurf(f.faces,f.vertices(:,1)/1000,f.vertices(:,2)/1000,f.vertices(:,3)/1000,...
            'edgecolor','none','facecolor',[57 110 154]/255,'facealpha',1,...
                'ambientstrength',light_ambient,...
                'diffusestrength',light_diffuse,...
                'specularstrength',light_specular,...
                'specularexponent',light_specularexponent);
        set(gca,'clim',[0.001 0.01])
        
%== OBJECTS FOR SCALE =====================================================


%     t2=trisurf(f_SCALE2.faces,...
%         f_SCALE2.vertices(:,1)/1000+2,f_SCALE2.vertices(:,2)/1000-8,f_SCALE2.vertices(:,3)/1000,...
%         'edgecolor','none','facecolor',[180 197 204]/255,'facealpha',1,...
%             'ambientstrength',light_ambient,...
%             'diffusestrength',light_diffuse,...
%             'specularstrength',light_specular,...
%             'specularexponent',light_specularexponent);
%         rotate(t2,[0 1 0],180)
% 
%     t3=trisurf(f_SCALE.faces,...
%         f_SCALE.vertices(:,1)/1000+2,f_SCALE.vertices(:,2)/1000+8,f_SCALE.vertices(:,3)/1000,...
%         'edgecolor','none','facecolor',[255 79 0]/255,'facealpha',1,...
%             'ambientstrength',light_ambient,...
%             'diffusestrength',light_diffuse,...
%             'specularstrength',light_specular,...
%             'specularexponent',light_specularexponent);
%         rotate(t3,[0 1 0],180)
%         
    t2=trisurf(f_SCALE2.faces,...
        f_SCALE2.vertices(:,1)/1000,f_SCALE2.vertices(:,2)/1000,f_SCALE2.vertices(:,3)/1000,...
        'edgecolor','none','facecolor',[180 197 204]/255,'facealpha',1,...
            'ambientstrength',light_ambient,...
            'diffusestrength',light_diffuse,...
            'specularstrength',light_specular,...
            'specularexponent',light_specularexponent);
        
        xxx=t2.Vertices;
        rotate(t2,[0 1 0],180,[mean(xxx(3,:)) mean(xxx(3,:)) mean(xxx(3,:))])
        xxx=t2.Vertices;
        rotate(t2,[0 0 1],180,[mean(xxx(3,:)) mean(xxx(3,:)) mean(xxx(3,:))])
        
        xxx=t2.Vertices;
        xxx(:,1)=xxx(:,1)-22;
        xxx(:,2)=xxx(:,2)-3;
        set(t2,'Vertices',xxx);

    t3=trisurf(f_SCALE.faces,...
        f_SCALE.vertices(:,1)/1000,f_SCALE.vertices(:,2)/1000,f_SCALE.vertices(:,3)/1000,...
        'edgecolor','none','facecolor',[255 79 0]/255,'facealpha',1,...
            'ambientstrength',light_ambient,...
            'diffusestrength',light_diffuse,...
            'specularstrength',light_specular,...
            'specularexponent',light_specularexponent);
        
        xxx=t3.Vertices;
        rotate(t3,[0 1 0],180,[mean(xxx(3,:)) mean(xxx(3,:)) mean(xxx(3,:))])
        xxx=t3.Vertices;
        rotate(t3,[0 0 1],90,[mean(xxx(3,:)) mean(xxx(3,:)) mean(xxx(3,:))])
        
        xxx=t3.Vertices;
        xxx(:,1)=xxx(:,1)-22;
        set(t3,'Vertices',xxx);
        
    %-- other things
    
        xxx=[-22 -6 0];
        scl=2;
        
        yyy=xxx;
        yyy(1)=yyy(1)+scl;
        mArrow3(xxx,yyy,'stemWidth',0.1,'tipWidth',0.2,'color','w',...
            'ambientstrength',light_ambient,...
            'diffusestrength',light_diffuse,...
            'specularstrength',light_specular,...
            'specularexponent',light_specularexponent);

        yyy=xxx;
        yyy(2)=yyy(2)+scl;
        mArrow3(xxx,yyy,'stemWidth',0.1,'tipWidth',0.2,'color','w',...
            'ambientstrength',light_ambient,...
            'diffusestrength',light_diffuse,...
            'specularstrength',light_specular,...
            'specularexponent',light_specularexponent);
        
        yyy=xxx;
        yyy(3)=yyy(3)+scl;
        mArrow3(xxx,yyy,'stemWidth',0.1,'tipWidth',0.2,'color','w',...
            'ambientstrength',light_ambient,...
            'diffusestrength',light_diffuse,...
            'specularstrength',light_specular,...
            'specularexponent',light_specularexponent);

        yyy=xxx;
        [xl,yl,zl]=sph2cart(125.2294*d2r,-61.8582*d2r,scl); % CA06      
        yyy(1)=yyy(1)+xl;     
        yyy(2)=yyy(2)+yl;    
        yyy(3)=yyy(3)+zl;
        mm1=mArrow3(xxx,yyy,'stemWidth',0.1,'tipWidth',0.2,'color','y',...
            'ambientstrength',0.5,...
            'diffusestrength',light_diffuse,...
            'specularstrength',light_specular,...
            'specularexponent',light_specularexponent,...
            'backfacelighting','unlit');
        
        plotcube([1 1 1],[-22 -10 0],1,[255 255 255]/255);
        

    %-- axis options
        daspect([1 1 1]);
        view([0 -90])
        axis vis3d;
        
        box on;
        grid on;
        spc=1;
        set(gca,'xtick',-30:spc:30,'ytick',-30:spc:30,'ztick',-30:spc:30)
        set(gca,'xticklabel',{},'yticklabel',{},'zticklabel',{});
        set(gca,'ticklen',[0 0])
        
        
exportfilename='global';
for view_option=0
    cl=camlight;
    delete(cl);

    if view_option==0
    % CA06 ----------------------------------------------------------------
        view([-176.7042+90, -51.0114]);
        camup([0 0 -1])
        camva(6);
        camroll(70);
        [xl,yl,zl]=sph2cart(125.2294*d2r,-61.8582*d2r,100); % CA06            
        cl=light('position',[xl yl zl],'style','infinite');
        camproj('ortho');
        axis off;
        VIEW_LABEL='CA06';
        set(gcf,'color','k');
        
    elseif view_option==1
    % VENTRAL -------------------------------------------------------------
        view([0 -89.999]);
        camva(6);
        camup([0 0 -1]);
        cl=camlight;
        VIEW_LABEL='VENTRAL';
        
        
    elseif view_option==2
    % DORSAL --------------------------------------------------------------
        view([0 89.999]);
        camva(6);
        camup([0 0 -1]);
        cl=camlight;
        VIEW_LABEL='DORSAL';

        
    elseif view_option==3
    % STARBOARD -----------------------------------------------------------
        view([0 0]);
        camva(6);
        camup([0 0 -1]);
        cl=camlight;
        VIEW_LABEL='STARBOARD';

        
    elseif view_option==4
    % PORT ----------------------------------------------------------------
        view([180 0]);
        camva(6);
        camup([0 0 -1]);
        cl=camlight;
        VIEW_LABEL='PORT';

        
    elseif view_option==5
    % FORE ----------------------------------------------------------------
        view([90 0]);
        camva(6);
        camup([0 0 -1]);
        cl=camlight;
        VIEW_LABEL='FORE';

        
    elseif view_option==6
    % AFT -----------------------------------------------------------------
        view([-90 0]);
        camva(6);
        camup([0 0 -1]);
        cl=camlight;
        VIEW_LABEL='AFT';
    end

axis tight;
    
    dosave=0 ;
    if dosave==1
        EXPORTFILENAME=['MU69_201211_',exportfilename,'_QUICKVIEW2_',VIEW_LABEL,'_SCALE'];
        EXPORTRESOLUTION='-r600';
        print(EXPORTFILENAME,'-dtiff',EXPORTRESOLUTION)
    end
end
%% SLICING THE NECK

figure(1);clf;hold on;
set(gcf,'color','w');
fonts=10;
fontn='helvetica neue';
set(gca,'fontsize',fonts,'fontname',fontn);

cc=[207 233 198;
    247 207 220]/255;
colormap(cc);

    %-- lighting properties (gentle)
        light_ambient=0.55;           % fractional strength of ambient light
        light_diffuse=0.4;           % fractional strength of diffuse light
        light_specular=0.05;          % fractional strength of specular light
        light_specularexponent=1;    % specular light exponent
        
        colormap(parula(100));
        [~,curvature] = findPointNormals(centers,[],[],true);

    %-- trisurf
        t1=trisurf(f.faces,f.vertices(:,1)/1000,f.vertices(:,2)/1000,f.vertices(:,3)/1000,...
            'edgecolor','none','cdata',curvature,'facealpha',0.5,...
            'facecolor',[1 1 1]*0.8,...
                'ambientstrength',light_ambient,...
                'diffusestrength',light_diffuse,...
                'specularstrength',light_specular,...
                'specularexponent',light_specularexponent);
        set(gca,'clim',[0.001 0.01])
    
%== FITTING NECK SLICE
    
    %-- planar slice
        neckfittype=fittype(@(a,b,c,y,z) a*y + b*z + c,...
                   'coefficients',{'a','b','c'},...
                   'independent',{'y','z'},...
                   'dependent','x');
        

%== PRESCRIBING NECK SLICE: MERGED MODEL, THINNEST CROSS-SECTION
    neckfit.a =    -0.13;
    neckfit.b =     0.06;
    neckfit.c =     -4.18;
        
    %-- plotting fit
        ytest=-8:0.1:8;
        [ytest,ztest]=meshgrid(ytest,ytest);
        xtest=neckfit.a*ytest + neckfit.b*ztest + neckfit.c;
        
    IN_NECK=inpolyhedron(f,[xtest(:) ytest(:) ztest(:)]*1000);
    
    FUCKITY=xtest.*0;
    FUCKITY(IN_NECK)=1;
  
    
        neck_surf=surf(xtest,ytest,ztest,'facealpha',0.5,...
            'edgecolor','none','cdata',FUCKITY);
        
        cc=colormap(parula(3));

    %-- axis options
        daspect([1 1 1]);
        %view([0 -90])
        view([-30 -30]);
        camup([0 0 -1]);
        axis vis3d;
        
        box on;
        grid on;
        spc=1;
        set(gca,'xtick',-30:spc:30,'ytick',-30:spc:30,'ztick',-30:spc:30)
        set(gca,'ticklen',[0 0])
        %xlabel('x (km)');ylabel('y (km)');zlabel('z (km)')
        set(gca,'xticklabel',{},'yticklabel',{},'zticklabel',{});
        
        %cb=colorbar('location','eastoutside');
        %set(cb,'fontsize',fonts,'fontname',fontn);

        camlight;
        
        pause(1);
    
%% CREATING HIGH-RESOLUTION MASCON GRID FOR CALCULATING VOLUME AND MOMENTS
disp(' ');disp('creating mascon grid for calculating moments, volumes, etc.');

%-- checking to see if the high-resolution grid has already been calculated
%     check = isfile(['MOMENTSGRID_201030_',global_mesh,'.mat']);
    check = isfile(['MOMENTSGRID_191126_',global_mesh,'.mat']);

if check==1
    disp(' (previous iteration found and loaded)')
    load(['MOMENTSGRID_191126_',global_mesh,'.mat'])
    yspc=xspc;
    zspc=xspc;
else

tic;
%-- making discrete grid
    xspc=100; %
    yspc=xspc;
    zspc=xspc;

%-- shape file vertices
    x=f.vertices(:,1);
    y=f.vertices(:,2);
    z=f.vertices(:,3);

%-- linear arrays, spanning the shape file
    xx=-max(-x(:)):xspc:max(x(:));
    yy=-max(-y(:)):yspc:max(y(:));
    zz=-max(-z(:)):zspc:max(z(:));

%-- discrete test point array, filling the entire volume
    [XX,YY,ZZ]=meshgrid(xx,yy,zz); 
    testpts=[XX(:), YY(:), ZZ(:)];
    IN=inpolyhedron(f,testpts);

%-- saving results
    save(['MOMENTSGRID_201030_',global_mesh,'.mat'],'IN','xspc','XX','YY','ZZ','testpts');
    
toc
end

%% CALCULATING MOMENTS OF INERTIA, CENTERS OF MASS, ETC.
disp(' ');disp('calculating moments of inertia, volumes, etc.');

MASS_ARRAY=zeros(3,1);
VOL_ARRAY =zeros(3,1);
COM_ARRAY =zeros(3,3);
MAX_PRINCIPAL_AXIS_ARRAY=zeros(3,3);
INT_PRINCIPAL_AXIS_ARRAY=zeros(3,3);
MIN_PRINCIPAL_AXIS_ARRAY=zeros(3,3);

figure(1);clf;hold on;
set(gcf,'color','w');
fonts=10;
fontn='helvetica neue';
set(gca,'fontsize',fonts,'fontname',fontn);
    
    %-- lighting properties
        light_ambient=0.1;           % fractional strength of ambient light
        light_diffuse=1.0;           % fractional strength of diffuse light
        light_specular=0.1;          % fractional strength of specular light
        light_specularexponent=1;    % specular light exponent
        
    [~,curvature] = findPointNormals(centers,[],[],true);

    %-- trisurf
        t1=trisurf(f.faces,f.vertices(:,1)/1000,f.vertices(:,2)/1000,f.vertices(:,3)/1000,...
            'edgecolor','none','facecolor','w','facealpha',0.5,...
                'ambientstrength',light_ambient,...
                'diffusestrength',light_diffuse,...
                'specularstrength',light_specular,...
                'specularexponent',light_specularexponent);
        set(gca,'clim',[0.001 0.01])
        camlight;   
    
    %-- axis options
        daspect([1 1 1]);
        view([0 -90])
        axis vis3d;
        
        box on;
        grid on;
        spc=1;
        set(gca,'xtick',-30:spc:30,'ytick',-30:spc:30,'ztick',-30:spc:30)
        set(gca,'ticklen',[0 0])
        xlabel('x (km)');ylabel('y (km)');zlabel('z (km)')
        

for BODY=[1 2 3]
% for BODY=3

    disp(' ');
    if BODY==1
        disp(' BIG LOBE:');
    elseif BODY==2
        disp(' SMALL LOBE:');
    elseif BODY==3
        disp(' ENTIRE BODY:');
    end
    
    if BODY==1
            in=IN;
            
            xtest=XX/1000;
            ytest=YY/1000;
            ztest=ZZ/1000;
            xtest_plane = (neckfit.a*ytest + neckfit.b*ztest + neckfit.c);
            LOBE_PTS=zeros(numel(IN),1);
            LOBE_PTS(xtest_plane<=xtest)=1;
            
            in1=in+LOBE_PTS;
            in=in1==2;
            
    elseif BODY==2
            in=IN;
            
            xtest=XX/1000;
            ytest=YY/1000;
            ztest=ZZ/1000;
            xtest_plane = (neckfit.a*ytest + neckfit.b*ztest + neckfit.c);
            LOBE_PTS=zeros(numel(IN),1);
            LOBE_PTS(xtest_plane>xtest)=1;
            
            in1=in+LOBE_PTS;
            in=in1==2;
            
    elseif BODY==3
            in=IN;
    end
    
    %-- mass array
        density=235;
        masses=double(in)*xspc*yspc*zspc*density;      % prisms
        volumes=double(in)*xspc*yspc*zspc;             % prisms
        
        MASS_ARRAY(BODY,1)=sum(masses);
        VOL_ARRAY(BODY,1)=sum(volumes);
        
    %-- storing the high resolution model volume for later
        high_resolution_volume = sum(volumes);
       
        disp([' volume = ',num2str(sum(volumes)/1000^3,'%5.4e'),' km^3'])
        disp([' volume equivalent diameter = ',num2str( (3/4/pi*sum(volumes))^(1/3)/1000*2,'%4.3f'),' km'])
        disp([' mass = ',num2str(sum(masses),'%4.3e'),' kg, assuming density of ',num2str(density),' kg/m^3']);
            
        disp([' x: max = ',num2str(max(XX(in))/1000,'%3.2f'),' km, min = ',num2str(min(XX(in))/1000,'%3.2f'),' km, delta = ',num2str((max(XX(in))-min(XX(in)))/1000,'%3.2f'),' km']);
        disp([' y: max = ',num2str(max(YY(in))/1000,'%3.2f'),' km, min = ',num2str(min(YY(in))/1000,'%3.2f'),' km, delta = ',num2str((max(YY(in))-min(YY(in)))/1000,'%3.2f'),' km']);
        disp([' z: max = ',num2str(max(ZZ(in))/1000,'%3.2f'),' km, min = ',num2str(min(ZZ(in))/1000,'%3.2f'),' km, delta = ',num2str((max(ZZ(in))-min(ZZ(in)))/1000,'%3.2f'),' km']);

       
    %-- center of mass
        COM=zeros(3,1);
        for i=1:1:numel(masses)
            COM=COM+masses(i)*[XX(i);YY(i);ZZ(i)];
        end
        COM=COM/sum(masses);
       
       COM_ARRAY(BODY,:)=COM';
       
       disp([' COM = [',num2str(COM(1)/1000,'%4.3f'),', ',num2str(COM(2)/1000,'%4.3f'),', ',num2str(COM(3)/1000,'%4.3f'),'] km']);
       
    %-- moments of inertia
        IXX=sum(masses(in).*((testpts(in,2)-COM(2)).^2+(testpts(in,3)-COM(3)).^2));
        IYY=sum(masses(in).*((testpts(in,1)-COM(1)).^2+(testpts(in,3)-COM(3)).^2));
        IZZ=sum(masses(in).*((testpts(in,1)-COM(1)).^2+(testpts(in,2)-COM(2)).^2));
       
        IXY=-sum(masses(in).*(testpts(in,1)-COM(1)).*(testpts(in,2)-COM(2)));
        IXZ=-sum(masses(in).*(testpts(in,1)-COM(1)).*(testpts(in,3)-COM(3)));
        IYZ=-sum(masses(in).*(testpts(in,2)-COM(2)).*(testpts(in,3)-COM(3)));
        
        I = [IXX IXY IXZ;
             IXY IYY IYZ;
             IXZ IYZ IZZ];
         
         period=15.9*60*60;
         omega=[0;0;2*pi/period];
         K = 1/2*omega'*I*omega;
         
    %-- determining principal axes
        [I_VECTORS,I_VALUES] = eig(I);
        I_VALUES=sum(I_VALUES);

    %-- determining maximum, intermediate, minimum eigenvalues and eigenvectors
        IND = [1;2;3];
        I_MAX = IND(I_VALUES==max(I_VALUES));
        I_MIN = IND(I_VALUES==min(I_VALUES));
            IND_INT=IND;
            IND_INT(I_MIN)=0;
            IND_INT(I_MAX)=0;
        I_INT=IND(sum(IND_INT));

        MAX_EIGENVECTOR = I_VECTORS(:,I_MAX);
        MIN_EIGENVECTOR = I_VECTORS(:,I_MIN);
        INT_EIGENVECTOR = cross(MAX_EIGENVECTOR,MIN_EIGENVECTOR);
        
        MAX_PRINCIPAL_AXIS_ARRAY(BODY,:)=MAX_EIGENVECTOR;
        INT_PRINCIPAL_AXIS_ARRAY(BODY,:)=INT_EIGENVECTOR;
        MIN_PRINCIPAL_AXIS_ARRAY(BODY,:)=MIN_EIGENVECTOR;
        
        disp([' moments of inertia: A = ',num2str(I_VALUES(I_MIN),'%3.2e'),...
              ', B = ',num2str(I_VALUES(I_INT),'%3.2e'),...
              ', C = ',num2str(I_VALUES(I_MAX),'%3.2e'),...
              ' kg m^2']);
          
       disp([' A=[',num2str(MIN_EIGENVECTOR(1),'%6.6f'),';',num2str(MIN_EIGENVECTOR(2),'%6.6f'),';',num2str(MIN_EIGENVECTOR(3),'%6.6f'),']'])
       disp([' B=[',num2str(INT_EIGENVECTOR(1),'%6.6f'),';',num2str(INT_EIGENVECTOR(2),'%6.6f'),';',num2str(INT_EIGENVECTOR(3),'%6.6f'),']'])
       disp([' C=[',num2str(MAX_EIGENVECTOR(1),'%6.6f'),';',num2str(MAX_EIGENVECTOR(2),'%6.6f'),';',num2str(MAX_EIGENVECTOR(3),'%6.6f'),']'])
       
end

%% POLYHEDRA GRAVITY CALCULATION
disp(' ');disp('performing polyhedra gravity calculation');

%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%
% Calculating polyhedron gravity after Robert A. Werner
% "The gravitational potential of a homogeneous polyhedron or don't cut
% corners," Celestial Mechanics and Dynamical Astronomy, 50, 253-279, 1994,
% http://adsabs.harvard.edu/abs/1994CeMDA..59..253W
%
% Original code developed by:
% Jim Richardson, Lunar and Planetary Laboratory, University of Arizona
% Version 1.0: 9 January 2003 ('gravmap')
%
% Modernized and modularized by:
% David Minton, Purdue University
% Version 2.0: 5 June 2007 ('gravmod.f90')
%
% Adapted for Python using f2py by:
% David Minton, Purdue University
% Version 3.0: 5 June 2019
%
% Adapted for MATLAB by:
% James Tuttle Keane, California Institute of Technology
% Version 4.0: 22 November 2019
%
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

%-- small number, used when quantities hit singularities
    smallnumber=1e-10;
    
%== SETTING UP TEST POINT GRID

    % GLOBAL SURFACE MESH - - - - - - - - - - - - - - - - - - - - - - - - -
    % (do this if you want to evaluate things on the surface of the body)
    
        testpoints=centers;
        
    % EXTERNAL GRID - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % (do this if you want to evaluate things in space around the body)
    
    
%         spc=0.2;
%         test_x = (-35:spc:35)*1000;
%         test_y = (-35:spc:35)*1000;
%         test_z = (-15:spc:15)*1000;
%         [TEST_X,TEST_Y,TEST_Z]=meshgrid(test_x,test_y,test_z);
%         testpoints=[TEST_X(:) TEST_Y(:) TEST_Z(:)];

        % NOTE: to reform the matrix, use reshape:
        %TEST_XR=reshape(testpoints(:,1),size(TEST_X,1),size(TEST_X,2),size(TEST_X,3));

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        

%== CREATING BLANK ARRAYS
    m=size(testpoints,1);
    gravitational_potential = zeros(m,1);
    gx = zeros(m,1);
    gy = zeros(m,1);
    gz = zeros(m,1);

%== LOOPING OVER POLYGON FACES
t_start=tic;
printoutat=round(linspace(1,n,20));
printoutat=printoutat(2:end-1);
printoutat=round(printoutat);
disp('  beginning loop')
for i=1:1:n
    %disp(i/n)
    
    %-- retrieving polygon vertices:
        x1 = vert1(i,1);
        y1 = vert1(i,2);
        z1 = vert1(i,3);
        x2 = vert2(i,1);
        y2 = vert2(i,2);
        z2 = vert2(i,3);
        x3 = vert3(i,1);
        y3 = vert3(i,2);
        z3 = vert3(i,3);
    
    %-- retrieving polygon normal:
        nhatx = normals(i,1);
        nhaty = normals(i,2);
        nhatz = normals(i,3);
    
    %-- calculate vector between vertices 1 and 2:
        x12 = x2 - x1;
        y12 = y2 - y1;
        z12 = z2 - z1;

    %-- calculate vector between vertices 2 and 3:
        x23 = x3 - x2;
        y23 = y3 - y2;
        z23 = z3 - z2;
    
    %-- defining the polygon reference frame, ihat = x12/|x12|:
        x12mag = sqrt(x12^2 + y12^2 + z12^2);
        x12mag(x12mag==0)=smallnumber;
        ihatx = x12 / x12mag;
        ihaty = y12 / x12mag;
        ihatz = z12 / x12mag;
    
    %-- defining the polygon reference frame, jhat = nhat x ihat:
        jhatx = (nhaty * ihatz) - (nhatz * ihaty);
        jhaty = (nhatz * ihatx) - (nhatx * ihatz);
        jhatz = (nhatx * ihaty) - (nhaty * ihatx);
        jmag = sqrt(jhatx^2 + jhaty^2 + jhatz^2);
        jmag(jmag==0)=smallnumber;
        jhatx = jhatx / jmag;
        jhaty = jhaty / jmag;
        jhatz = jhatz / jmag;

    %-- calculating the positions of vertices in the polygon reference 
    %   frame (the components of the three basis vectors, ihat, jhat, nhat, 
    %   form the rotation matrix to transform from body-centric reference 
    %   frame to the polygon reference frame):
        xi1   = (ihatx*x1) + (ihaty*y1) + (ihatz*z1);
        eta1  = (jhatx*x1) + (jhaty*y1) + (jhatz*z1);
        zeta1 = (nhatx*x1) + (nhaty*y1) + (nhatz*z1);
        xi2   = (ihatx*x2) + (ihaty*y2) + (ihatz*z2);
        eta2  = (jhatx*x2) + (jhaty*y2) + (jhatz*z2);
        zeta2 = (nhatx*x2) + (nhaty*y2) + (nhatz*z2);
        xi3   = (ihatx*x3) + (ihaty*y3) + (ihatz*z3);
        eta3  = (jhatx*x3) + (jhaty*y3) + (jhatz*z3);
        zeta3 = (nhatx*x3) + (nhaty*y3) + (nhatz*z3);
    
    %-- calculating the lengths of each polygon edge:
    
        %-  side 1-2:
            r12 = sqrt((xi2 - xi1)^2 + (eta2 - eta1)^2 + (zeta2 - zeta1)^2);
            r12(r12==0)=smallnumber;
            
        %-  side 2-3:
            r23 = sqrt((xi3 - xi2)^2 + (eta3 - eta2)^2 + (zeta3 - zeta2)^2);
            r23(r23==0)=smallnumber;
            
        %-  side 3-1:
            r31 = sqrt((xi1 - xi3)^2 + (eta1 - eta3)^2 + (zeta1 - zeta3)^2);
            r31(r31==0)=smallnumber;
        
    %-- rotating test points to the polygon reference frame:
        x0p = (ihatx*testpoints(:,1)) + (ihaty*testpoints(:,2)) + (ihatz*testpoints(:,3));
        y0p = (jhatx*testpoints(:,1)) + (jhaty*testpoints(:,2)) + (jhatz*testpoints(:,3));
        z0p = (nhatx*testpoints(:,1)) + (nhaty*testpoints(:,2)) + (nhatz*testpoints(:,3));
    
    %-- calculate difference between test points and polygon vertices:
        dx1 = xi1   - x0p;
        dy1 = eta1  - y0p;
        dz1 = zeta1 - z0p;
        dx2 = xi2   - x0p;
        dy2 = eta2  - y0p;
        dz2 = zeta2 - z0p;
        dx3 = xi3   - x0p;
        dy3 = eta3  - y0p;
        dz3 = zeta3 - z0p;
        dz = (dz1 + dz2 + dz3)/3;
        dz(dz==0)=smallnumber;
         
    %-- distance from test point to vertices:
        r1 = ((dx1.*dx1)+(dy1.*dy1)+(dz1.*dz1)).^(1/2);
        r1(r1==0)=smallnumber;
        
        r2 = ((dx2.*dx2)+(dy2.*dy2)+(dz2.*dz2)).^(1/2);
        r2(r2==0)=smallnumber;
        
        r3 = ((dx3.*dx3)+(dy3.*dy3)+(dz3.*dz3)).^(1/2);
        r3(r3==0)=smallnumber;
        
    %-- calculating determinant and logarithm expressions (equations on the
    %   top of p.263 of Werner 1994):
        det12 = (dx1 .* dy2) - (dx2 .* dy1);
        det23 = (dx2 .* dy3) - (dx3 .* dy2);
        det31 = (dx3 .* dy1) - (dx1 .* dy3);
        
        % note: 'd' = denominator, 'n' = numerator.
        
        dl = (r1 + r2 - r12);
        dl(dl==0)=smallnumber;
        nl = (r1 + r2 + r12) ./ dl;
        nl(nl<=0)=smallnumber;
        L12 = log(nl);
        L12 = L12 ./ r12;
        
        dl = (r2 + r3 - r23);
        dl(dl==0)=smallnumber;
        nl = (r2 + r3 + r23) ./ dl;
        nl(nl<=0)=smallnumber;
        L23 = log(nl);
        L23 = L23 ./ r23;
        
        dl = (r3 + r1 - r31);
        dl(dl==0)=smallnumber;
        nl = (r3 + r1 + r31) ./ dl;
        nl(nl<=0)=smallnumber;
        L31 = log(nl);
        L31 = L31 ./ r31;
        
    %-- calculating the components in the arctangent terms (the 'S' terms 
    %   on p.263 of Werner 1994):
    
        % note: 'denom' = denominator, 'num' = numerator.
        
        nnum = (xi1*(eta2-eta3))+(xi2*(eta3-eta1))+(xi3*(eta1-eta2));
        denom1 = ((xi2-xi1)*(xi1-xi3))+((eta2-eta1)*(eta1-eta3));
        denom2 = ((xi3-xi2)*(xi2-xi1))+((eta3-eta2)*(eta2-eta1));
        denom3 = ((xi1-xi3)*(xi3-xi2))+((eta1-eta3)*(eta3-eta2));
         
    %-- outer brackets and S terms on p.263
         snum = dz * nnum;
         sdenom1 = ((det31.*det12) + ((dz.^2)*denom1)) ./ (-1.0d0*r1);
         sdenom1(sdenom1==0)=smallnumber;
         
         sdenom2 = ((det12.*det23) + ((dz.^2)*denom2)) ./ (-1.0d0*r2);
         sdenom2(sdenom2==0)=smallnumber;
         
         sdenom3 = ((det23.*det31) + ((dz.^2)*denom3)) ./ (-1.0d0*r3);
         sdenom3(sdenom3==0)=smallnumber;
        
         S1 = atan(snum ./ sdenom1);
             rng1=snum>0;
             rng2=sdenom1<0;
             rng3=rng1+rng2;
             rng=rng3==2;
             S1(rng)=pi+S1(rng);
             
             rng1=snum<0;
             rng2=sdenom1<0;
             rng3=rng1+rng2;
             rng=rng3==2;
             S1(rng)=S1(rng)-pi;
         
         S2 = atan(snum ./ sdenom2);
             rng1=snum>0;
             rng2=sdenom2<0;
             rng3=rng1+rng2;
             rng=rng3==2;
             S2(rng)=pi+S2(rng);
             
             rng1=snum<0;
             rng2=sdenom2<0;
             rng3=rng1+rng2;
             rng=rng3==2;
             S2(rng)=S2(rng)-pi;
             
         S3 = atan(snum ./ sdenom3);
             rng1=snum>0;
             rng2=sdenom3<0;
             rng3=rng1+rng2;
             rng=rng3==2;
             S3(rng)=pi+S3(rng);
             
             rng1=snum<0;
             rng2=sdenom3<0;
             rng3=rng1+rng2;
             rng=rng3==2;
             S3(rng)=S3(rng)-pi;
             
    %-- sum factor with sign-dependent pi
        sfac=dz.*0;
        sfac(dz>=0)= 1.0d0;
        sfac(dz<0) =-1.0d0;
        af = (S1 + S2 + S3) - (sfac*pi);
        
    %-- calculating the gravitational potential from the polyhedron:
        dU = (dz/2) .* (det12.*L12 + det23.*L23 + det31.*L31) - (dz.^2/2 .* af);
        
        % note, the test points cannot be inside the body, else the
        % derivation breaks down
         
    %-- calculating the gravitational acceleration, in the polygon 
    %   reference frame (the 'S' terms on p.264-265 of Werner 1994):
        dUdx0 = -dz/2 .* ((eta2-eta1)*L12 + (eta3-eta2)*L23 + (eta1-eta3)*L31);
        
        dUdy0 =  dz/2 .* ((xi2-xi1)*L12 + (xi3-xi2)*L23 + (xi1-xi3)*L31);
        
        dUdz0 = -1/2  .* (det12.*L12 + det23.*L23 + det31.*L31) + dz.*af;
        
    %-- converting the accelerations from the polygon reference frame to 
    %   the body-centric reference frame (note that inverse rotation matrix
    %   is just transpose of rotation matrix):
        dUdx = (ihatx*dUdx0) + (jhatx*dUdy0) + (nhatx*dUdz0);
        dUdy = (ihaty*dUdx0) + (jhaty*dUdy0) + (nhaty*dUdz0);
        dUdz = (ihatz*dUdx0) + (jhatz*dUdy0) + (nhatz*dUdz0);
         
    %-- finally, summing the total potentials and accelerations:
        gravitational_potential = gravitational_potential + dU;
        gx  = gx  + dUdx;
        gy  = gy  + dUdy;
        gz  = gz  + dUdz;
    
        
    t_run=toc(t_start);
    t_step=t_run/i;
    t_remaining=(n-i)*t_step;
    if sum(i==printoutat)>=1
    disp(['  t-minus: ',num2str(t_remaining/60,'%3.2f'),' minutes'])
    end
end

%-- printing run time
    t_run=toc(t_start);
    disp(['  total run time: ',num2str(t_run/60,'%3.2f'),' minutes'])
    disp('  calculations complete');
    
    
%-- saving the output (so you don't have to do this step again)
    save('MU69_201201_GravityFrom_GlobalModel_EvaluatedOn_GlobalModel.mat',...
        'gx','gy','gz','dUdx','dUdy','dUdz','gravitational_potential',...
        'testpoints');

%% IMPORTING PREVIOUSLY SOLVED SOLUTION

    load('MU69_201201_GravityFrom_MergedModel_EvaluatedOn_MergedModel')
    
%% CALCULATING GRAVITY AND ROTATIONAL POTENTIALS AND ACCELERATIONS
disp(' ');disp('calculating gravity+rotational potentials and accelerations')


center_of_mass  = COM_ARRAY(3,:);
rotation_vector = MAX_PRINCIPAL_AXIS_ARRAY(3,:);

rotation_period = 15.92*60*60;           % seconds (MU69)
rotation_rate   = 2*pi/rotation_period;  % radians/second

rotation_vector=rotation_vector/norm(rotation_vector);

%-- calculating the distance from the test points to the rotation vector:

    rrr = ((center_of_mass(1)-testpoints(:,1))*rotation_vector(1) + ...
    (center_of_mass(2)-testpoints(:,2))*rotation_vector(2) + ...
    (center_of_mass(3)-testpoints(:,3))*rotation_vector(3));

    rp_x = (center_of_mass(1)-testpoints(:,1)) - ...
           rrr*rotation_vector(1);
    rp_y = (center_of_mass(2)-testpoints(:,2)) - ...
           rrr*rotation_vector(2);
    rp_z = (center_of_mass(3)-testpoints(:,3)) - ...
           rrr*rotation_vector(3);

    rp = (rp_x.^2 + rp_y.^2 + rp_z.^2).^(1/2);
    
%-- calculating centrifugal potential and acceleration
    ROTATIONAL_POTENTIAL = 1/2 * rotation_rate^2 * rp.^2;
    
    CX = rotation_rate^2*-rp_x;
    CY = rotation_rate^2*-rp_y;
    CZ = rotation_rate^2*-rp_z;
    
    CC = (CX.^2 + CY.^2 + CZ.^2).^(1/2);

%-- scaling the polyhedron model by the density of the object:
    
    density=235; % kg/m^3
%     density=155; % kg/m^3
%     density=500; % kg/m^3
%     density=2000; % kg/m^3
   
    GRAVITATIONAL_POTENTIAL = gravitational_potential * G * density;
    
    GX = gx * G * density;
    GY = gy * G * density;
    GZ = gz * G * density;
    
    GG = (GX.^2 + GY.^2 + GZ.^2).^(1/2);

%-- combining rotational and gravitational components:
    AX = GX + CX;
    AY = GY + CY;
    AZ = GZ + CZ;
    
    
    TOTAL_ACCELERATION_MAGNITUDE=(AX.^2 + AY.^2 + AZ.^2).^(1/2);
    
%% CALCULATING QUANTITIES ON THE SURFACE OF THE MESH
disp(' ');disp('calculating surface forces, potentials, slopes')

%-- calculating radial accelerations:
    A_RADIAL_X = normals(:,1) .* AX;
    A_RADIAL_Y = normals(:,2) .* AY;
    A_RADIAL_Z = normals(:,3) .* AZ;
    
    A_RADIAL = (A_RADIAL_X.^2 + A_RADIAL_Y.^2 + A_RADIAL_Z.^2).^(1/2);
    
FORCE_RADIAL=zeros(n,1);
FORCE_MAG=zeros(n,1);
SLOPE=zeros(n,1);
NORTH=zeros(n,3);
EAST=zeros(n,3);
FORCE_EAST=zeros(n,1);
FORCE_NORTH=zeros(n,1);
DIP_DIRECTION=zeros(n,1);
CORIOLIS=zeros(n,1);

for i=1:1:n
    
    %-- normals
        nv=normals(i,:);
        nv=nv/norm(nv);

    %-- radial force
        fv=[AX(i) AY(i) AZ(i)];
        FORCE_RADIAL(i) = dot(fv,nv);
        FORCE_MAG(i)=norm(fv);

    %-- surface slope:
        SLOPE(i) = 180-acosd(dot(fv/norm(fv),nv));
        
    %-- slope vector
%         SLOPE_VECTOR = fv/norm(fv)+nv/norm(nv);
%                 fv=[AX(i) AY(i) AZ(i)];
        
%         SLOPE_VECTOR = fv-dot(fv,nv);
%         
%         OMEGA = rotation_rate*rotation_vector;
%         
%         CORIOLIS(i)=-1*norm(cross(OMEGA,SLOPE_VECTOR));

    %-- local geodetic reference frame
        east=cross([0 0 1],nv);
        east=east/norm(east);
        north=cross(nv,east);
        north=north/norm(north);

        NORTH(i,:)=north;
        EAST(i,:)=east;

        FORCE_EAST(i)=dot(fv/norm(fv),east);
        FORCE_NORTH(i)=dot(fv/norm(fv),north);
        DIP_DIRECTION(i)=90-atan2d(FORCE_NORTH(i),FORCE_EAST(i));
    
end

%== GEOPOTENTIAL
    GEOPOTENTIAL=-ROTATIONAL_POTENTIAL-GRAVITATIONAL_POTENTIAL;

%== GEOMETRIC ALTITUDE (distance from the COM)
    GEOMETRIC_ALTITUDE = ((testpoints(:,1)-COM(1)).^2 + ...
                          (testpoints(:,2)-COM(2)).^2 + ...
                          (testpoints(:,3)-COM(3)).^2).^(1/2);
                  
%== GEOPOTENTIAL ALTITUDE
    REFERENCE_POTENTIAL=min(GEOPOTENTIAL(:));
    GEOPOTENTIAL_ALTITUDE = (GEOPOTENTIAL-REFERENCE_POTENTIAL)./TOTAL_ACCELERATION_MAGNITUDE;
    
%== JACOBI SPEED
%   Eq. 13 of Scheeres et al. 2016
%   Jacobi speed is the speed a particle would gain if it went from the
%   geopotential high to the geopotential low.
    JACOBI_SPEED = (-2*GEOPOTENTIAL).^(1/2);
    JACOBI_SPEED = (JACOBI_SPEED.^2 - min(JACOBI_SPEED).^2).^(1/2);
    
%== SURFACE ESCAPE SPEEDS
%   Eq. 10.13 (p.241) of Scheeres 2012
    ESCAPE_SPEED=zeros(n,1);
    for i=1:1:n
        
    %-- normals
        nv=normals(i,:);
%         nv=normals_STEREO(i,:);
        nv=nv/norm(nv);
        
    %-- radial force
        fv=[AX(i) AY(i) AZ(i)];
        
        ftv = fv-dot(nv,fv);
        ftv = ftv/norm(ftv);
        
    %-- vector from rotation vector to test point
        rp_vector = -[rp_x(i); rp_y(i); rp_z(i)];
        
    %-- maximum potential
        UMAX=GRAVITATIONAL_POTENTIAL(i);
        
    %-- rotation vector (scaled by rotation rate)
        rotation_vector2=rotation_vector*rotation_rate;
        
    %-- escape speed
        ESCAPE_SPEED(i) = -dot(ftv',cross(rotation_vector2',rp_vector)) + ...
                          (dot(ftv',cross(rotation_vector2',rp_vector))^2 + ...
                           2*UMAX - ...
                           dot(cross(rotation_vector2',rp_vector),cross(rotation_vector2',rp_vector))).^(1/2);
%         ESCAPE_SPEED(i) = -dot(nv',cross(rotation_vector2',rp_vector)) + ...
%                           (dot(nv',cross(rotation_vector2',rp_vector))^2 + ...
%                            2*UMAX - ...
%                            dot(cross(rotation_vector2',rp_vector),cross(rotation_vector2',rp_vector))).^(1/2);
    end

%== ROCHE LOBE
%   This is the energy of the minimum equilibrium point (and must be
%   calculated beyond the surface of the object)
    minimum_equilibrium_point_potential=  -5.714133737995875; %MU69 235 kg/m^3);
    
%== RETURN SPEED 
    RETURN_SPEED = (2*-(GEOPOTENTIAL-minimum_equilibrium_point_potential)).^(1/2);
    

%% PLOTTING RESULTS

for MEGA=9

hfig=figure(10);clf;hold on;
set(gcf,'color','w');
fonts=6;
fontn='helvetica neue';
set(gca,'fontsize',fonts,'fontname',fontn);

%-- figure size
    set(gcf,'Units','centimeters')    
    set(gcf,'position',[15 24 40 22.5]);
    set(gcf,'PaperUnits','centimeters');
    P0=get(gcf,'position');
    set(gcf,'PaperPosition',[0 0 P0(3) P0(4);])
    set(gcf,'inverthardcopy','off');
        
    %-- lighting properties (gentle)
        light_ambient=0.55;           % fractional strength of ambient light
        light_diffuse=0.4;           % fractional strength of diffuse light
        light_specular=0.05;          % fractional strength of specular light
        light_specularexponent=1;    % specular light exponent
        
%== DATASET ===============================================================

DO_SLOPE_VECTORS=0;

% Note, I use a lot of custom colormaps (e.g., 'wiec.txt'). You can replace
% all of those "importdata" lines with generic MATLAB colormaps (e.g.,
% cc=colormap(parula).
    
if MEGA==1
%-- WHITE
    cc=[1 1 1];
    colormap(cc)
    set(gca,'clim',[0 1]);
    Q=SLOPE;
    exportfilename='White';
    cb_label=' ';

elseif MEGA==2
%-- SLOPE
%     cc=importdata('coolwarm2.txt')/255;
    cc=importdata('wiec.txt')/255;
    colormap(cc)
    set(gca,'clim',[0 30]);
    Q=SLOPE;
    exportfilename='Slope';
    cb_label='slope, degrees';
    DO_SLOPE_VECTORS=0;

elseif MEGA==3
%-- ACCELERATIONS
    cc=importdata('wiec.txt')/255;
    colormap(cc);
    set(gca,'clim',[0 0.55]);
    Q=A_RADIAL*1000;
    exportfilename='RadialAcceleration';
    cb_label='total radial acceleration, mm/s^2';

elseif MEGA==4
%-- ACCELERATIONS
    cc=importdata('wiec.txt')/255;
    colormap(cc);
    set(gca,'clim',[0 0.55]);
    Q=TOTAL_ACCELERATION_MAGNITUDE*1000;
    exportfilename='TotalAcceleration';
    cb_label='total acceleration, mm/s^2';
    
%     Q=CC*1000;
%     exportfilename='CentrifugalAcceleration';
%     cb_label='centrifugal acceleration, mm/s^2';

%     Q=GG*1000;
%     exportfilename='GravitationalAcceleration';
%     cb_label='gravitational acceleration, mm/s^2';
    
%     Q=A_RADIAL*1000;
%     exportfilename='RadialAcceleration';
%     cb_label='total radial acceleration, mm/s^2';
    

elseif MEGA==5
%-- GEOPOTENTIAL
    cc=importdata('birefring3.txt')/255;
    cc=flipud(cc);
    colormap(cc);
    Q=GEOPOTENTIAL;
    set(gca,'clim',[min(Q(:)) max(Q(:))])
    set(gca,'clim',minimum_equilibrium_point_potential+[-0.53 0.53])
    exportfilename='Geopotential';
    cb_label='geopotential, m^2/s^2';

elseif MEGA==6
%-- GEOPOTENTIAL ALTITUDE
    cc=importdata('batlow.txt')/255;
    cc=brighten(cc,0.25);
    colormap(cc);
    Q=GEOPOTENTIAL_ALTITUDE/1000;
    set(gca,'clim',[min(Q(:)) max(Q(:))])
    exportfilename='GeopotentialAltitude2';
    cb_label='geopotential altitude, km';
    
elseif MEGA==7 
%-- GEOMETRIC ALTITUDE
    cc=importdata('batlow.txt')/255;
    colormap(cc);
    Q=GEOMETRIC_ALTITUDE/1000;
    set(gca,'clim',[min(Q(:)) max(Q(:))])
    exportfilename='GeometricAltitude';
    cb_label='geometric altitude, km';

elseif MEGA==8
%-- JACOBI SPEED
    cc=importdata('wiec.txt')/255;
    colormap(cc);
    Q=JACOBI_SPEED;
    set(gca,'clim',[0 1.5])
    exportfilename='JacobiSpeed';
    cb_label='Jacobi speed, m/s';

elseif MEGA==9
%-- ESCAPE SPEED
    cc=importdata('wiec.txt')/255;
    colormap(cc);
    Q=ESCAPE_SPEED;
    set(gca,'clim',[0 4.5])
    exportfilename='EscapeSpeed';
    cb_label='direct escape speed, m/s';

elseif MEGA==10
%-- RETURN SPEED
    cc=importdata('wiec.txt')/255;
    colormap(cc);
    Q=real(RETURN_SPEED);
    set(gca,'clim',[min(Q(:)) max(Q(:))])
    set(gca,'clim',[0 1.5])
    exportfilename='ReturnSpeed';
    cb_label='return speed, m/s';
    
elseif MEGA==11
%-- TEST
    cc=importdata('wiec.txt')/255;
    colormap(cc);
    Q=real(JACOBI_SPEED)./real(ESCAPE_SPEED);
    exportfilename='Test';
    cb_label='test';
    
end


%== COLORBAR ==============================================================

    cb=colorbar('location','eastoutside');
    set(cb,'tickdir','both','fontsize',fonts,'fontname',fontn);
    ylabel(cb,cb_label,'fontsize',fonts,'fontname',fontn);
    cb_pos=get(cb,'position');
    set(cb,'position',[0.92 cb_pos(2) cb_pos(3) cb_pos(4)]);
    
    set(cb,'color','k');
        

%== SURFACE PLOT ==========================================================

    trisurf(f.faces,f.vertices(:,1)/1000,f.vertices(:,2)/1000,f.vertices(:,3)/1000,...
        'edgecolor','none','cdata',Q,'facealpha',1,...
        'ambientstrength',light_ambient,...
        'diffusestrength',light_diffuse,...
        'specularstrength',light_specular,...
        'specularexponent',light_specularexponent);
%     
%     trisurf(f.faces,f.vertices(:,1)/1000,f.vertices(:,2)/1000,f.vertices(:,3)/1000,...
%         'edgecolor','none','facecolor','w','facealpha',1,...
%         'ambientstrength',light_ambient,...
%         'diffusestrength',light_diffuse,...
%         'specularstrength',light_specular,...
%         'specularexponent',light_specularexponent);
    
%== AXIS OPTIONS ==========================================================
        daspect([1 1 1]);
        box on;
        grid on;
        spc=1;
        set(gca,'xtick',-30:spc:30,'ytick',-30:spc:30,'ztick',-30:spc:30)
        set(gca,'xticklabel',{},'yticklabel',{},'zticklabel',{})
        set(gca,'ticklen',[0 0])
        %xlabel('x (km)');ylabel('y (km)');zlabel('z (km)')
        axis tight;
        axis vis3d;
        %axis off;
        
        
        camtarget(COM_ARRAY(3,:)/1000);
        
%== SPIN VECTOR ===========================================================

    %-- spin vector
        vec=COM_ARRAY(3,:)/1000+[0 0 -5.5];
        scl=1;
        clr=[134 6 37]/255;
        stemwidth=0.2;
        tipwidth=0.0;
        mArrow3(COM_ARRAY(3,:)/1000,...
                [vec(1)*scl vec(2)*scl vec(3)*scl],...
                'stemWidth',stemwidth,'tipWidth',tipwidth,...
                'color',clr,...
                'ambientstrength',light_ambient,...
                'diffusestrength',light_diffuse,...
                'specularstrength',light_specular,...
                'specularexponent',light_specularexponent);
            
            theta=linspace(0,2*pi,300);
            p1=patch(stemwidth*cos(theta)+vec(1),stemwidth*sin(theta)+vec(2),theta.*0+vec(3),'r');
            set(p1,'EdgeColor','none',...
                'ambientstrength',light_ambient,... 
                'diffusestrength',light_diffuse,...
                'specularstrength',light_specular,...
                'specularexponent',light_specularexponent);
        vec=COM_ARRAY(3,:)/1000+[0 0 6.5];
        scl=1;
        clr=[134 6 37]/255;
        stemwidth=0.2;
        tipwidth=0.3;
        mArrow3(COM_ARRAY(3,:)/1000,...
                [vec(1)*scl vec(2)*scl vec(3)*scl],...
                'stemWidth',stemwidth,'tipWidth',tipwidth,...
                'color',clr,...
                'ambientstrength',light_ambient,...
                'diffusestrength',light_diffuse,...
                'specularstrength',light_specular,...
                'specularexponent',light_specularexponent);
            

        
%== OBJECTS FOR SCALE =====================================================

    t2=trisurf(f_SCALE2.faces,...
        f_SCALE2.vertices(:,1)/1000+2,f_SCALE2.vertices(:,2)/1000-8,f_SCALE2.vertices(:,3)/1000,...
        'edgecolor','none','facecolor',[180 197 204]/255,'facealpha',1,...
            'ambientstrength',light_ambient,...
            'diffusestrength',light_diffuse,...
            'specularstrength',light_specular,...
            'specularexponent',light_specularexponent);
        rotate(t2,[0 1 0],180)

    t3=trisurf(f_SCALE.faces,...
        f_SCALE.vertices(:,1)/1000+2,f_SCALE.vertices(:,2)/1000+8,f_SCALE.vertices(:,3)/1000,...
        'edgecolor','none','facecolor',[255 79 0]/255,'facealpha',1,...
            'ambientstrength',light_ambient,...
            'diffusestrength',light_diffuse,...
            'specularstrength',light_specular,...
            'specularexponent',light_specularexponent);
        rotate(t3,[0 1 0],180)
            
%== DOWNSLOPE ARROWS ======================================================

if DO_SLOPE_VECTORS==1
    random_i=randperm(size(centers,1));
    lwidth=0.5;
    
%for i=random_i(1:1:4000)

% - - - - - - - - - - - - - - - - - (Z)
for x_test=[-20:0.5:17]*1000
    for y_test=[-11:0.5:11]*1000
        for fuck_test=[-1 1]
        [intersect,~,~,~,~]=...
        TriangleRayIntersection([x_test, y_test, 0*1000],...
                                [x_test, y_test, fuck_test*100*1000],...
                                vert1,vert2,vert3);
% % - - - - - - - - - - - - - - - - - (Y)
% for x_test=[-20:0.5:17]*1000
%     for z_test=[-5:0.5:5]*1000
%         for fuck_test=[-1 1]
%         [intersect,~,~,~,~]=...
%         TriangleRayIntersection([x_test, 0*1000, z_test],...
%                                 [x_test, fuck_test*100*1000, z_test],...
%                                 vert1,vert2,vert3);
% % - - - - - - - - - - - - - - - - - (X)
% for y_test=[-11:0.5:11]*1000
%     for z_test=[-5:0.5:5]*1000
%         for fuck_test=[1]
%         [intersect,~,~,~,~]=...
%         TriangleRayIntersection([10*1000, y_test z_test],...
%                                 [-1*100*1000, y_test z_test],...
%                                 vert1,vert2,vert3);
        sum(intersect);
        if sum(intersect)>=1
            
            facet_indices=1:1:n;
            ii=facet_indices(intersect);
            
            
        for i=ii


        north=NORTH(i,:);
        east =EAST(i,:);

        fuck=[FORCE_EAST(i); FORCE_NORTH(i)];
        fuck=fuck/norm(fuck);


        %scl_super=1;
        scl_super=SLOPE(i)/30*1.2;
        scl_super=max([0.4 scl_super]);

        VECTOR=east*FORCE_EAST(i)+north*FORCE_NORTH(i);
        VECTOR=VECTOR/norm(VECTOR) * 1.0*scl_super;
        plot3([centers(i,1)/1000 centers(i,1)/1000+VECTOR(1)],...
              [centers(i,2)/1000 centers(i,2)/1000+VECTOR(2)],...
              [centers(i,3)/1000 centers(i,3)/1000+VECTOR(3)],...
              '-k',...
              'linewidth',lwidth)

        rv=normals(i,:);
        t=-60*d2r;
        RM = [cos(-(t-pi/2)) + rv(1)^2*(1-cos(-(t-pi/2))),              rv(1)*rv(2)*(1-cos(-(t-pi/2))) - rv(3)*sin(-(t-pi/2)),      rv(1)*rv(3)*(1-cos(-(t-pi/2))) + rv(2)*sin(-(t-pi/2));...
             rv(2)*rv(1)*(1-cos(-(t-pi/2))) + rv(3)*sin(-(t-pi/2)),     cos(-(t-pi/2)) + rv(2)^2*(1-cos(-(t-pi/2))),                rv(2)*rv(3)*(1-cos(-(t-pi/2))) - rv(1)*sin(-(t-pi/2));...
             rv(3)*rv(1)*(1-cos(-(t-pi/2))) - rv(2)*sin(-(t-pi/2)),     rv(3)*rv(2)*(1-cos(-(t-pi/2))) + rv(1)*sin(-(t-pi/2)),      cos(-(t-pi/2)) + rv(3)^2*(1-cos(-(t-pi/2)))];
        vector=RM*VECTOR';
        vector=vector/norm(vector)*0.5*scl_super;
        plot3([centers(i,1)/1000+VECTOR(1) centers(i,1)/1000+VECTOR(1)+vector(1)],...
              [centers(i,2)/1000+VECTOR(2) centers(i,2)/1000+VECTOR(2)+vector(2)],...
              [centers(i,3)/1000+VECTOR(3) centers(i,3)/1000+VECTOR(3)+vector(3)],...
              '-k',...
              'linewidth',lwidth)

        rv=normals(i,:);
        t=240*d2r;
        RM = [cos(-(t-pi/2)) + rv(1)^2*(1-cos(-(t-pi/2))),              rv(1)*rv(2)*(1-cos(-(t-pi/2))) - rv(3)*sin(-(t-pi/2)),      rv(1)*rv(3)*(1-cos(-(t-pi/2))) + rv(2)*sin(-(t-pi/2));...
             rv(2)*rv(1)*(1-cos(-(t-pi/2))) + rv(3)*sin(-(t-pi/2)),     cos(-(t-pi/2)) + rv(2)^2*(1-cos(-(t-pi/2))),                rv(2)*rv(3)*(1-cos(-(t-pi/2))) - rv(1)*sin(-(t-pi/2));...
             rv(3)*rv(1)*(1-cos(-(t-pi/2))) - rv(2)*sin(-(t-pi/2)),     rv(3)*rv(2)*(1-cos(-(t-pi/2))) + rv(1)*sin(-(t-pi/2)),      cos(-(t-pi/2)) + rv(3)^2*(1-cos(-(t-pi/2)))];
        vector=RM*VECTOR';
        vector=vector/norm(vector)*0.5*scl_super;
        plot3([centers(i,1)/1000+VECTOR(1) centers(i,1)/1000+VECTOR(1)+vector(1)],...
              [centers(i,2)/1000+VECTOR(2) centers(i,2)/1000+VECTOR(2)+vector(2)],...
              [centers(i,3)/1000+VECTOR(3) centers(i,3)/1000+VECTOR(3)+vector(3)],...
              '-k',...
              'linewidth',lwidth)
        end
        end
        end
    end
end
end
    
%== CREATING IMAGES =======================================================


for view_option=[1:1:6]
delete(cl);

    if view_option==0
    % CA06 ----------------------------------------------------------------
        view([-176.7042+90, -51.0114]);
        camup([0 0 -1])
        camva(6);
        camroll(70);
        [xl,yl,zl]=sph2cart(125.2294*d2r,-61.8582*d2r,100); % CA06            
        light('position',[xl yl zl],'style','infinite')
        camproj('ortho');
        axis off;
        VIEW_LABEL='CA06';
        set(gcf,'color','k');
        axis off;
        set(cb,'color','w');
        
    elseif view_option==1
    % VENTRAL -------------------------------------------------------------
        view([0 -89.999]);
        camva(6);
        camup([0 0 -1]);
        cl=camlight;
        VIEW_LABEL='VENTRAL';
        set(gcf,'color','w');
        axis on;
        
        
    elseif view_option==2
    % DORSAL --------------------------------------------------------------
        view([0 89.999]);
        camva(6);
        camup([0 0 -1]);
        cl=camlight;
        VIEW_LABEL='DORSAL';
        set(gcf,'color','w');
        axis on;

        
    elseif view_option==3
    % STARBOARD -----------------------------------------------------------
        view([0 0]);
        camva(6);
        camup([0 0 -1]);
        cl=camlight;
        VIEW_LABEL='STARBOARD';
        set(gcf,'color','w');
        axis on;

        
    elseif view_option==4
    % PORT ----------------------------------------------------------------
        view([180 0]);
        camva(6);
        camup([0 0 -1]);
        cl=camlight;
        VIEW_LABEL='PORT';
        set(gcf,'color','w');
        axis on;

        
    elseif view_option==5
    % FORE ----------------------------------------------------------------
        view([90 0]);
        camva(6);
        camup([0 0 -1]);
        cl=camlight;
        VIEW_LABEL='FORE';
        set(gcf,'color','w');
        axis on;

        
    elseif view_option==6
    % AFT -----------------------------------------------------------------
        view([-90 0]);
        camva(6);
        camup([0 0 -1]);
        cl=camlight;
        VIEW_LABEL='AFT';
        set(gcf,'color','w');
        axis on;
        
    elseif view_option==100
        
    % WTF --------------------------------------------------------------
        view([-40 -27]);
        camva(6.5);
        camup([0 0 -1]);
        cl=camlight;
        VIEW_LABEL='WTF';
        set(gcf,'color','w');
        axis on;
    end
        
        xlim([-20 17])
        zlim([-6 7]);
        ylim([-11 11])
%         axis tight;
    
    dosave=0;
    if dosave==1
        EXPORTFILENAME=['MU69_201211_',exportfilename,'_',num2str(density,'%4.1i'),'_',VIEW_LABEL,'_wtf2'];
        EXPORTRESOLUTION='-r600';
        print(EXPORTFILENAME,'-dtiff',EXPORTRESOLUTION)
    end
end

%== MAKING ANIMATION ======================================================


do_vid=0;
    if do_vid==1
    spc=0.75;
    vid_counter=1;
    for theta=360:-spc:0+spc
        disp(['video: ',num2str(theta/(360-spc)*100),'% complete']);

%         [xc,yc,zc]=sph2cart((theta)*d2r,-40*d2r,1);
%         [xl,yl,zl]=sph2cart((theta+50)*d2r,-40*d2r,1);
        [xc,yc,zc]=sph2cart((theta)*d2r,-35*d2r,1);
        [xl,yl,zl]=sph2cart((theta+50)*d2r,-40*d2r,1);
        
        campos([xc yc zc]*250);
        camup([0 0 -1]);
        camtarget(COM/1000);
        camva(5.2);
        axis vis3d;
        axis off;
        vid_filename=['MU69_200930_DPS_',exportfilename];
        
        
        if vid_counter==1
            l1=light('position',[xl yl zl],'style','infinite');
                v=VideoWriter(vid_filename,'MPEG-4'); %#ok<TNMLP>
                v.Quality=95;
                open(v);
                frame=getframe(hfig);
                writeVideo(v,frame);

        else
            set(l1,'position',[xl yl zl]);
                frame=getframe(hfig);
                writeVideo(v,frame);
        end
        vid_counter=vid_counter+1;
    end
    close(v);
    else
    end
    
    
end
    