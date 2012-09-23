/*
 * panorama implementation
 *
 * Copyright (C) HamidReza Houshiar
 *
 * Released under the GPL version 3.
 *
 */

#include "slam6d/fbr/panorama.h"

using namespace std;

namespace fbr{
  
  void panorama::init(unsigned int width, unsigned int height, projection_method method, unsigned int numberOfImages, double param, panorama_map_method mMethod){
    iWidth = width;
    iHeight = height;
    pMethod = method;
    nImages = numberOfImages;
    pParam = param;
    if(mMethod == FARTHEST){
      iMap.create(iHeight, iWidth, CV_32FC(3));
      iMap = cv::Scalar::all(0);
    }
    else if(mMethod == EXTENDED){
      extendedIMap.resize(iHeight);
      for (unsigned int i = 0; i < iHeight; i++)
	extendedIMap[i].resize(iWidth);
    }
    iReflectance.create(iHeight, iWidth, CV_8U);
    iReflectance = cv::Scalar::all(0);
    iRange.create(iHeight, iWidth, CV_32FC(1));
    iRange = cv::Scalar::all(0);
    mapMethod = mMethod;
  }

  panorama::panorama(unsigned int width, unsigned int height, projection_method method, unsigned int numberOfImages, double param, panorama_map_method mMethod){ 
    init(width, height, method, numberOfImages, param, mMethod);
  }
  
  panorama::panorama(unsigned int width, unsigned int height, projection_method method, unsigned int numberOfImages, double param){ 
    init(width, height, method, numberOfImages, param, FARTHEST);
  }
  
  panorama::panorama(unsigned int width, unsigned int height, projection_method method, unsigned int numberOfImages){
    double param = 0;
    if(method == PANNINI)
      param = 1;
    else if (method == STEREOGRAPHIC)
      param = 2;
   
    init(width, height, method, numberOfImages, param, FARTHEST);
  }  

  panorama::panorama(unsigned int width, unsigned int height, projection_method method){ 
    double param = 0;
    unsigned int numberOfImages = 1;
    if(method == RECTILINEAR)
      numberOfImages = 3;
    else if(method == PANNINI){
      numberOfImages = 3;
      param = 1;
    } else if (method == STEREOGRAPHIC){
      numberOfImages = 3;
      param = 2;
    }

    init(width, height, method, numberOfImages, param, FARTHEST);
  }
  
  void panorama::map(int x, int y, cv::MatIterator_<cv::Vec4f> it, double range){
    iReflectance.at<uchar>(y,x) = (*it)[3]*255;//reflectance
    iRange.at<float>(y,x) = range;//range
    if(mapMethod == FARTHEST){
      //adding the point with max distance
      if( iRange.at<float>(y,x) < range ){
	iMap.at<cv::Vec3f>(y,x)[0] = (*it)[0];//x
	iMap.at<cv::Vec3f>(y,x)[1] = (*it)[1];//y
	iMap.at<cv::Vec3f>(y,x)[2] = (*it)[2];//z
      }
    }else if(mapMethod == EXTENDED){
      //adding all the points
      cv::Vec3f point;
      point[0] = (*it)[0];//x
      point[1] = (*it)[1];//y
      point[2] = (*it)[2];//z
      extendedIMap[y][x].push_back(point);
    }
  }

  void panorama::createPanorama(cv::Mat scan){

    //EQUIRECTANGULAR projection
    if(pMethod == EQUIRECTANGULAR){
      //adding the longitude to x axis and latitude to y axis
      double xFactor = (double) iWidth / 2 / M_PI;
      int widthMax = iWidth - 1;
      double yFactor = (double) iHeight / ((MAX_ANGLE - MIN_ANGLE) / 360 * 2 * M_PI);
      //shift all the valuse to positive points on image 
      double heightLow =(0 - MIN_ANGLE) / 360 * 2 * M_PI;
      int heightMax = iHeight - 1;
      
      cv::MatIterator_<cv::Vec4f> it, end; 
      
      for( it = scan.begin<cv::Vec4f>(), end = scan.end<cv::Vec4f>(); it != end; ++it){
	double kart[3], polar[3], phi, theta, range;
	kart[0] = (*it)[2]/100;
	kart[1] = (*it)[0]/-100;
	kart[2] = (*it)[1]/100;
	toPolar(kart, polar);
	//theta == polar[0] == scan [4]
	//phi == polar[1] == scan [5]
	//range == polar[2] == scan [3]
	theta = polar[0] * 180 / M_PI;
	phi = polar[1] * 180 / M_PI;
	range = polar[2];
	//horizantal angle of view of [0:360] and vertical of [-40:60]
	phi = 360.0 - phi;
	phi = phi * 2.0 * M_PI / 360.0;
	theta -= 90;
	theta *= -1;
	theta *= 2.0 * M_PI / 360.0;
	int x = (int) ( xFactor * phi);
	if (x < 0) x = 0;
	if (x > widthMax) x = widthMax;
	int y = (int) ( yFactor * (theta + heightLow) );
	y = heightMax - y;
	if (y < 0) y = 0;
	if (y > heightMax) y = heightMax;
	
	//create the iReflectance iRange and map
	map(x, y, it, range);
      }
    }
    
    //CYLINDRICAL projection
    if(pMethod == CYLINDRICAL){
      //adding the longitude to x and tan(latitude) to y
      //find the x and y range
      double xFactor = (double) iWidth / 2 / M_PI;
      int widthMax = iWidth - 1;
      double yFactor = (double) iHeight / (tan(MAX_ANGLE / 360 * 2 * M_PI) - tan(MIN_ANGLE / 360 * 2 * M_PI));
      double heightLow = (MIN_ANGLE) / 360 * 2 * M_PI;
      int heightMax = iHeight - 1;
      
      cv::MatIterator_<cv::Vec4f> it, end; 
      
      for( it = scan.begin<cv::Vec4f>(), end = scan.end<cv::Vec4f>(); it != end; ++it){
	double kart[3], polar[3], phi, theta, range;
	kart[0] = (*it)[2]/100;
	kart[1] = (*it)[0]/-100;
	kart[2] = (*it)[1]/100;
	toPolar(kart, polar);
	//theta == polar[0] == scan [4]
	//phi == polar[1] == scan [5]
	//range == polar[2] == scan [3]
	theta = polar[0] * 180 / M_PI;
	phi = polar[1] * 180 / M_PI;
	range = polar[2];
	//horizantal angle of view of [0:360] and vertical of [-40:60]
	phi = 360.0 - phi;
	phi = phi * 2.0 * M_PI / 360.0;
	theta -= 90;
	theta *= -1;
	theta *= 2.0 * M_PI / 360.0;
	int x = (int) ( xFactor * phi);
	if (x < 0) x = 0;
	if (x > widthMax) x = widthMax;
	int y = (int) ((double) yFactor * (tan(theta) - tan(heightLow)));
	y = heightMax - y;
	if (y < 0) y = 0;
	if (y > heightMax) y = heightMax;
	
	//create the iReflectance iRange and map
	map(x, y, it, range);
      }
    }
    
    //Mercator Projection
    if( pMethod == MERCATOR){
      //find the x and y range
      double xFactor = (double) iWidth / 2 / M_PI;
      int widthMax = iWidth - 1;
      double yFactor = (double) iHeight / ( log( tan( MAX_ANGLE / 360 * 2 * M_PI ) + ( 1 / cos( MAX_ANGLE / 360 * 2 * M_PI ) ) ) - log ( tan( MIN_ANGLE / 360 * 2 * M_PI) + (1/cos(MIN_ANGLE / 360 * 2 * M_PI) ) ) );
      double heightLow = log(tan(MIN_ANGLE / 360 * 2 * M_PI) + (1/cos(MIN_ANGLE / 360 * 2 * M_PI)));
      int heightMax = iHeight - 1;
      
      cv::MatIterator_<cv::Vec4f> it, end; 
      
      for( it = scan.begin<cv::Vec4f>(), end = scan.end<cv::Vec4f>(); it != end; ++it){
	double kart[3], polar[3], phi, theta, range;
	kart[0] = (*it)[2]/100;
	kart[1] = (*it)[0]/-100;
	kart[2] = (*it)[1]/100;
	toPolar(kart, polar);
	//theta == polar[0] == scan [4]
	//phi == polar[1] == scan [5]
	//range == polar[2] == scan [3]
	theta = polar[0] * 180 / M_PI;
	phi = polar[1] * 180 / M_PI;
	range = polar[2];
	//horizantal angle of view of [0:360] and vertical of [-40:60]
	phi = 360.0 - phi;
	phi = phi * 2.0 * M_PI / 360.0;
	theta -= 90;
	theta *= -1;
	theta *= 2.0 * M_PI / 360.0;
	int x = (int) ( xFactor * phi);
	if (x < 0) x = 0;
	if (x > widthMax) x = widthMax;
	int y = (int) ( yFactor * (log(tan(theta) + (1/cos(theta))) - heightLow) );
	y = heightMax - y;
	if (y < 0) y = 0;
	if (y > heightMax) y = heightMax;
	
	//create the iReflectance iRange and map
	map(x, y, it, range);
      }
    }
    
    //RECTILINEAR projection
    if(pMethod == RECTILINEAR){
      //default value for nImages
      if(nImages == 0) nImages = 3;
      cout<<"Number of images per scan is: "<<nImages<<endl;
      double l0, p1, iMinx, iMaxx, iMiny, iMaxy, interval;
      interval = 2 * M_PI / nImages;
      iMiny = -M_PI/9;
      iMaxy = 2*M_PI/9;
      //latitude of projection center
      p1 = 0;
      
      //go through all points 
      cv::MatIterator_<cv::Vec4f> it, end; 
      
      for( it = scan.begin<cv::Vec4f>(), end = scan.end<cv::Vec4f>(); it != end; ++it){
	double kart[3], polar[3], phi, theta, range;
	kart[0] = (*it)[2]/100;
	kart[1] = (*it)[0]/-100;
	kart[2] = (*it)[1]/100;
	toPolar(kart, polar);
	//theta == polar[0] == scan [4]
	//phi == polar[1] == scan [5]
	//range == polar[2] == scan [3]
	theta = polar[0] * 180 / M_PI;
	phi = polar[1] * 180 / M_PI;
	range = polar[2];
	//horizantal angle of view of [0:360] and vertical of [-40:60]
	phi = 360.0 - phi;
	phi = phi * 2.0 * M_PI / 360.0;
	theta -= 90;
	theta *= -1;
	theta *= 2.0 * M_PI / 360.0;
	for(unsigned int j = 0 ; j < nImages ; j++){
	  iMinx = j * interval;
	  iMaxx = (j + 1) * interval;
	  //check for point in interval
	  if(phi < iMaxx && phi > iMinx){
	    double max, min, coscRectilinear;
	    //the longitude of projection center
	    l0 = iMinx + interval / 2;
	    //finding the min and max of the x direction
	    coscRectilinear = sin(p1) * sin(iMaxy) + cos(p1) * cos(iMaxy) * cos(iMaxx - l0);
	    max = (cos(iMaxy) * sin(iMaxx - l0) / coscRectilinear);
	    coscRectilinear = sin(p1) * sin(iMiny) + cos(p1) * cos(iMiny) * cos(iMinx - l0);
	    min = (cos(iMiny) * sin(iMinx - l0) / coscRectilinear);
	    double xFactor = (double) (iWidth / nImages) / (max - min);
	    double xlow = min;
	    int widthMax = (iWidth / nImages) - 1;
	    //finding the min and max of y direction
	    coscRectilinear = sin(p1) * sin(iMaxy) + cos(p1) * cos(iMaxy) * cos(iMaxx - l0);
	    max = ( (cos(p1) * sin(iMaxy) - sin(p1) * cos(iMaxy) * cos(iMaxx - l0) )/ coscRectilinear);
	    coscRectilinear = sin(p1) * sin(iMiny) + cos(p1) * cos(iMiny) * cos(iMinx - l0);
	    min = ( (cos(p1) * sin(iMiny) - sin(p1) * cos(iMiny) * cos(iMinx - l0) )/ coscRectilinear);
	    double yFactor = (double) iHeight / (max - min);
	    double heightLow = min;
	    int heightMax = iHeight - 1;
	    //project the points and add them to image
	    coscRectilinear = sin(p1) * sin(theta) + cos(p1) * cos(theta) * cos(phi - l0);
	    int x = (int)(xFactor) * ((cos(theta) * sin(phi - l0) / coscRectilinear) - xlow);
	    if (x < 0) x = 0;
	    if (x > widthMax) x = widthMax;
	    x = x + (j * iWidth / nImages);
	    int y = (int) (yFactor) * (( (cos(p1) * sin(theta) - sin(p1) * cos(theta) * cos(phi - l0)) / coscRectilinear) - heightLow);
	    y = heightMax - y;
	    if (y < 0) y = 0;
	    if (y > heightMax) y = heightMax;
	    
	    //create the iReflectance iRange and map
	    map(x, y, it, range);
	  }
	}
      }
    }
    
    //PANNINI projection
    if(pMethod == PANNINI){
      //default values for nImages and dPannini==pParam
      if(pParam == 0) pParam = 1;
      if(nImages == 0) nImages = 3;
      cout << "Parameter d is:" << pParam <<", Number of images per scan is:" << nImages << endl;
      double l0, p1, iMinx, iMaxx, iMiny, iMaxy, interval;
      interval = 2 * M_PI / nImages;
      iMiny = -M_PI/9;
      iMaxy = 2*M_PI/9;
      //latitude of projection center
      p1 = 0;
      
      cv::MatIterator_<cv::Vec4f> it, end; 
         
      for( it = scan.begin<cv::Vec4f>(), end = scan.end<cv::Vec4f>(); it != end; ++it){
	double kart[3], polar[3], phi, theta, range;
	kart[0] = (*it)[2]/100;
	kart[1] = (*it)[0]/-100;
	kart[2] = (*it)[1]/100;
	toPolar(kart, polar);
	//theta == polar[0] == scan [4]
	//phi == polar[1] == scan [5]
	//range == polar[2] == scan [3]
	theta = polar[0] * 180 / M_PI;
	phi = polar[1] * 180 / M_PI;
	range = polar[2];
	//horizantal angle of view of [0:360] and vertical of [-40:60]
	phi = 360.0 - phi;
	phi = phi * 2.0 * M_PI / 360.0;
	theta -= 90;
	theta *= -1;
	theta *= 2.0 * M_PI / 360.0;
	for(unsigned int j = 0 ; j < nImages ; j++){
	  iMinx = j * interval;
	  iMaxx = (j + 1) * interval;
	  //check for point in interval
	  if(phi < (iMaxx) && phi > (iMinx)){
	    double max, min, sPannini;
	    //the longitude of projection center
	    l0 = iMinx + interval / 2;
	    //use the S variable of pannini projection mentioned in the thesis
	    //finding the min and max of the x direction
	    sPannini = (pParam + 1) / (pParam + sin(p1) * tan(iMaxy) + cos(p1) * cos(iMaxx - l0));
	    max = sPannini * (sin(iMaxx - l0));
	    sPannini = (pParam + 1) / (pParam + sin(p1) * tan(iMiny) + cos(p1) * cos(iMinx - l0));
	    min = sPannini * (sin(iMinx - l0));
	    double xFactor = (double) (iWidth / nImages) / (max - min);
	    double xlow = min;
	    int widthMax = (iWidth / nImages) - 1;
	    //finding the min and max of y direction
	    sPannini = (pParam + 1) / (pParam + sin(p1) * tan(iMaxy) + cos(p1) * cos(iMaxx - l0));
	    max = sPannini * (tan(iMaxy) * (cos(p1) - sin(p1) * 1/tan(iMaxy) * cos(iMaxx - l0)));
	    sPannini = (pParam + 1) / (pParam + sin(p1) * tan(iMiny) + cos(p1) * cos(iMinx - l0));
	    min = sPannini * (tan(iMiny) * (cos(p1) - sin(p1) * 1/tan(iMiny) * cos(iMinx - l0)));
	    double yFactor = (double) iHeight / (max - min);
	    double heightLow = min;
	    int heightMax = iHeight - 1;
	    //project the points and add them to image
	    sPannini = (pParam + 1) / (pParam + sin(p1) * tan(theta) + cos(p1) * cos(phi - l0));
	    int x = (int)(xFactor) * (sPannini * sin(phi - l0) - xlow);
	    if (x < 0) x = 0;
	    if (x > widthMax) x = widthMax;
	    x = x + (j * iWidth / nImages);
	    int y = (int) (yFactor) * ( (sPannini * tan(theta) * (cos(p1) - sin(p1) * (1/tan(theta)) * cos(phi - l0) ) ) - heightLow );
	    y = heightMax - y;
	    if (y < 0) y = 0;
	    if (y > heightMax) y = heightMax;
		    
	    //create the iReflectance iRange and map
	    map(x, y, it, range);
	  }
	}
      }
    }
    
    //STEREOGRAPHIC projection
    if(pMethod == STEREOGRAPHIC){
      //default values for nImages and rStereographic==pParam
      if(pParam == 0) pParam = 2;
      if(nImages == 0) nImages = 3;
      cout << "Paremeter R is:" << pParam << ", Number of images per scan is:" << nImages << endl;
      // l0 and p1 are the center of projection iminx, imaxx, iminy, imaxy are the bounderis of intervals
      double l0, p1, iMinx, iMaxx, iMiny, iMaxy, interval;
      interval = 2 * M_PI / nImages;
      iMiny = -M_PI/9;
      iMaxy = 2*M_PI/9;
      //latitude of projection center
      p1 = 0;
      
      //go through all points
      cv::MatIterator_<cv::Vec4f> it, end; 
         
      for( it = scan.begin<cv::Vec4f>(), end = scan.end<cv::Vec4f>(); it != end; ++it){
	double kart[3], polar[3], phi, theta, range;
	kart[0] = (*it)[2]/100;
	kart[1] = (*it)[0]/-100;
	kart[2] = (*it)[1]/100;
	toPolar(kart, polar);
	//theta == polar[0] == scan [4]
	//phi == polar[1] == scan [5]
	//range == polar[2] == scan [3]
	theta = polar[0] * 180 / M_PI;
	phi = polar[1] * 180 / M_PI;
	range = polar[2];
	//horizantal angle of view of [0:360] and vertical of [-40:60]
	phi = 360.0 - phi;
	phi = phi * 2.0 * M_PI / 360.0;
	theta -= 90;
	theta *= -1;
	theta *= 2.0 * M_PI / 360.0;
	for (unsigned int j = 0 ; j < nImages ; j++){
	  iMinx = j * interval;
	  iMaxx = (j + 1) * interval;
	  //check for point in intervals
	  if(phi < (iMaxx) && phi > (iMinx)){
	    double max, min, k;
	    //longitude of projection center
	    l0 = iMinx + interval / 2;
	    //use the R variable of stereographic projection mentioned in the thesis
	    //finding the min and max of x direction
	    k = (2 * pParam) / (1 + sin(p1) * sin(p1) + cos(p1) * cos(p1) * cos(iMaxx - l0));
	    max = k * cos(p1) * sin (iMaxx - l0);
	    k = (2 * pParam) / (1 + sin (p1) * sin(p1) + cos(p1) * cos(p1) * cos(iMinx -l0));
	    min = k * cos(p1) * sin (iMinx -l0);
	    double xFactor = (double) (iWidth / nImages) / (max - min);
	    double xlow = min;
	    int widthMax = (iWidth / nImages) - 1;
	    //finding the min and max of y direction
	    k = (2 * pParam) / (1 + sin(p1) * sin(iMaxy) + cos(p1) * cos(iMaxy) * cos(iMaxx - l0));
	    max = k * (cos(p1) * sin(iMaxy) - sin(p1) * cos(iMaxy) * cos(iMaxx - l0));
	    k = (2 * pParam) / (1 + sin(p1) * sin(iMiny) + cos(p1) * cos(iMiny) * cos(iMinx - l0));
	    min = k * (cos(p1) * sin(iMiny) - sin(p1) * cos(iMiny) * cos(iMinx - l0));
	    double yFactor = (double) iHeight / (max - min);
	    double heightLow = min;
	    int heightMax = iHeight - 1;
	    //project the points and add them to image
	    k = (2 * pParam) / (1 + sin(p1) * sin(theta) + cos(p1) * cos(theta) * cos(phi - l0));
	    int x = (int) (xFactor) * (k * cos(theta) * sin(phi - l0) - xlow);
	    if (x < 0) x = 0;
	    if (x > widthMax) x = widthMax;
	    x = x + (j * iWidth / nImages);
	    int y = (int) (yFactor) * (k * ( cos(p1) * sin(theta) - sin(p1) * cos(theta) * cos(phi - l0) ) - heightLow);
	    y = heightMax - y;
	    if (y < 0) y = 0;
	    if (y > heightMax) y = heightMax;
	    
	    //create the iReflectance iRange and map
	    map(x, y, it, range);
	  }
	}
      }
    }
    
    //ZAXIS projection
    if(pMethod == ZAXIS){
      double zmin = -200;
      double zmax = 4000;
      //adding the longitude to x axis and latitude to y axis
      double xFactor = (double) iWidth / 2 / M_PI;
      int widthMax = iWidth - 1;
      cout << "ZMAX= " << zmax << " ZMIN= "<< zmin << endl;
      double yFactor = (double) iHeight / (zmax - zmin);
      //shift all the valuse to positive points on image 
      double heightLow = zmin;
      int heightMax = iHeight - 1;
      
      cv::MatIterator_<cv::Vec4f> it, end; 
         
      for( it = scan.begin<cv::Vec4f>(), end = scan.end<cv::Vec4f>(); it != end; ++it){
	double kart[3], polar[3], phi, theta, range;
	kart[0] = (*it)[2]/100;
	kart[1] = (*it)[0]/-100;
	kart[2] = (*it)[1]/100;
	toPolar(kart, polar);
	//theta == polar[0] == scan [4]
	//phi == polar[1] == scan [5]
	//range == polar[2] == scan [3]
	theta = polar[0] * 180 / M_PI;
	phi = polar[1] * 180 / M_PI;
	range = polar[2];
	//horizantal angle of view of [0:360] and vertical of [-40:60]
	phi = 360.0 - phi;
	phi = phi * 2.0 * M_PI / 360.0;
	theta -= 90;
	theta *= -1;
	theta *= 2.0 * M_PI / 360.0;
	int x = (int) ( xFactor * phi);
	if (x < 0) x = 0;
	if (x > widthMax) x = widthMax;
	///////////////////check this
	int y = (int) ( yFactor * ((*it)[1] - heightLow) );
	y = heightMax - y;
	if (y < 0) y = 0;
	if (y > heightMax) y = heightMax;
	
	//create the iReflectance iRange and map
	map(x, y, it, range);
      }
    }

    // nomenclature is chosen according to:
    // http://mathworld.wolfram.com/AlbersEqual-AreaConicProjection.html
    if(pMethod == CONIC){
      double phi, lambda, theta, rho;
      int widthMax = iWidth - 1;
      int heightMax = iHeight - 1;

      // the following code is not tested with phi_0 and lambda_0 of different
      // values than zero
      double phi_0 = 0.0;
      double lambda_0 = 0.0;

      // the panorama class does not allow for passing of arguments, so the
      // following values of phi_1 and phi_2 for the standard parallels are
      // fixed during runtime.
      double phi_1 = 0.0;
      double phi_2 = M_PI / 3.0; // 60 degrees
      // Nevertheless, the following code allows for phi_1 and phi_2 to be
      // different from the default values above

      // the angle over which the arc that the projection creates extends is
      // calculated with by
      //
      //  tau = (sin(phi_1)+sin(phi_2))*pi
      //
      // eg for phi_1 = 0.0 and phi_2 = pi/2 = 90 deg, the arc would be a half circle
      //    for phi_1 = 0.0 and phi_2 = pi/6 = 30 deg, the arc would be a quarter circle
      //
      // the variable n stores the value of tau/(2*pi)

      double n = 0.5*(sin(phi_1) + sin(phi_2));
      double C = cos(phi_1)*cos(phi_1) + 2 * n * sin(phi_1);
      double rho_0 = sqrt(C - 2 * n * sin(phi_0))/n;

      double scale, xshift, yshift, pwidth, pheight, xscale, yscale;
      // the value n decides whether the arc is less than, equal or more than
      // 180 degrees
      if (n >= 0.5) {
        // arc is bigger than or equal to 180 degrees

        // the x shift is equal to the radius which can be found by the
        // y coordinate of the projection of polar coordinate (-pi/2, 0.0)
        xshift = 2.0-rho_0+(sqrt(C+2*n)/n)*cos(-n*lambda_0);
        // the width is twice the radius
        pwidth = 2.0*xshift;

        // the y shift is equal to the y coordinate of the projection of polar
        // coordinate (-pi/2, pi)
        yshift = rho_0-(sqrt(C+2*n)/n)*cos(n*(M_PI-lambda_0));

        // the height is equal to the difference of the y coordinate of the
        // projection of polar coordinate (-pi/2, pi) and the y coordinate of
        // the projection of polar coordinate (-pi/2, 0.0)
        pheight = (sqrt(C+2*n)/n)*(cos(-n*lambda_0)-cos(n*(M_PI-lambda_0)));
      } else {
        // arc is smaller than 180 degrees

        // the x shift is equal to the x coordinate of the projection of polar
        // coordinate (-pi/2, pi)
        xshift = (sqrt(C+2*n)/n)*sin(n*(M_PI-lambda_0));
        // the width is twice the x shift
        pwidth = 2.0*xshift;

        // the y shift is equal to the y coordinate of the projection of polar
        // coordinate (pi/2, pi)
        yshift = rho_0-(sqrt(C-2*n)/n)*cos(n*(M_PI-lambda_0));

        // the height is equal to the difference of the y coordinate of the
        // projection of polar coordinate (pi/2, pi) and the y coordinate of
        // the projection of polar coordinate (-pi/2, 0.0)
        pheight = (sqrt(C+2*n)/n)*cos(-n*lambda_0)-(sqrt(C-2*n)/n)*cos(n*(M_PI-lambda_0));
      }

      xscale = iWidth / pwidth;
      yscale = iHeight / pheight;
      // pick the smaller scaling
      if (xscale < yscale)
          scale = xscale;
      else
          scale = yscale;

      cv::MatIterator_<cv::Vec4f> it, end;

      for( it = scan.begin<cv::Vec4f>(), end = scan.end<cv::Vec4f>(); it != end; ++it){
        double kart[3], polar[3], range;
        kart[0] = (*it)[2]/100;
        kart[1] = (*it)[0]/-100;
        kart[2] = (*it)[1]/100;
        toPolar(kart, polar);
        //theta == polar[0] == scan [4]
        //phi == polar[1] == scan [5]
        //range == polar[2] == scan [3]

        // considered turning phi upside down (negating it), so that the
        // "bottom" of the scan is on the longer side (bottom) of the cone
        // but the results do not look well as the scan is squeezed into a
        // small stripe at the bottom (at least with the default parallels)
        // and thus not well visible
        //
        // convert polar[0] from [0..pi] to [-pi/2..pi/2]
        // convert polar[1] from [0..2pi] to [-pi..pi]
        phi = polar[0] - M_PI/2.0;
        lambda = polar[1] - M_PI;
        range = polar[2];
        theta = n * (lambda - lambda_0);
        rho = sqrt(C - 2 * n * sin(phi))/n;
        // before shifting and scaling the center of the arc is at x,y = (0.0,
        // 1.0) and looks like this for n < 0.5:
        //
        //                       A Y
        //                       |
        //                       |
        //                       + (0,1)
        //       (pi/2,-pi) *    |    * (pi/2,pi)
        //                       * (pi/2,0)
        //(-pi/2,-pi) *          |          * (-pi/2,pi)
        //       ----------------+----------------->
        //                       |                  X
        //                       * (-pi/2,0)
        //                       |
        //
        int x = (int) ((rho * sin(theta) + xshift)*scale);
        // turn upside down because image y axis points downward not upward
        int y = (int) ((-1)*(rho_0 - rho * cos(theta) - yshift)*scale);
        if (x < 0) x = 0;
        if (x > widthMax) x = widthMax;
        if (y < 0) y = 0;
        if (y > heightMax) y = heightMax;
        map(x, y, it, range);
      }
    }
  }

  unsigned int panorama::getImageWidth(){
    return iWidth;
  }

  unsigned int panorama::getImageHeight(){
    return iHeight;
  }

  projection_method panorama::getProjectionMethod(){
    return pMethod;
  }

  unsigned int panorama::getNumberOfImages(){
    return nImages;
  }
  
  double panorama::getProjectionParam(){
    return pParam;
  }
  
  cv::Mat panorama::getReflectanceImage(){
    return iReflectance;
  }

  cv::Mat panorama::getMap(){
    return iMap;
  }

  cv::Mat panorama::getRangeImage(){
    return iRange;
  }
  
  vector<vector<vector<cv::Vec3f> > > panorama::getExtendedMap(){
    return extendedIMap;
  }

  panorama_map_method panorama::getMapMethod(){
    return mapMethod;
  }

  void panorama::getDescription(){
    cout<<"panorama created with width: "<<iWidth<<", and height: "<<iHeight<<", and projection method: "<<projectionMethodToString(pMethod)<<", number of images: "<<nImages<<", projection param: "<<pParam<<"."<<endl;
    cout<<endl;
  }
}

