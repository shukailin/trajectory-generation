#include <iostream>
#include <cmath>
#include <ros/ros.h>
#include <geometry_msgs/Twist.h>
#include <traj_generation/FlightCmd.h>

#define PI 3.14159
#define GRAVITY 9.8    //m/sec^2
#define MASS_OF_QUAD 1.5 //kg

double z[4] = {0.5, -0.5, 0.5, -0.5};
double x[4] = {0.5, -0.5, 0.5, -0.5};

double gLog[1000] = {0};

typedef struct
{
  double initPos[4]; //{x_0, y_0, z_0, yaw_0}
  double initVel[4]; //{Vx_0, Vy_0, Vz_0, Vyaw_0}
  double goalPos[4]; //{x_f, y_f, z_f, yaw_f}
  double goalVel[4]; //{Vx_f, Vy_f, Vz_f, Vyaw_f}
}STrajInitParams;


typedef struct
{
  double qt[4]; //{x, y, z, yaw}
  double qdt[4]; //{Vx, Vy, Vz, Vyaw}
  double qddt[4]; //{Ax, Ay, Az, Ayaw}
}STrajParams;

typedef struct
{
  double A[4];
  double B[4];
  double C[4];
  double E[4];
  double F[4];
}SPolyParams;

/*
  Output command Structure for Flight Control
*/
typedef struct
{
  float roll;
  float pitch;
  float yaw;
  float throttle;
}SFlightCmd;


class CTrajGeneration
{
  public:
	CTrajGeneration(ros::NodeHandle &n);        
	bool TestFunction(void);
	void TrajLoop( void );
  private:
	float m_samplingTime;	
	ros::Publisher pub;
	geometry_msgs::Twist baseCmd;
  traj_generation::FlightCmd m_flightCmd;
	SPolyParams sPolyParams;
  STrajParams sTrajParams;
  SFlightCmd sFlightCmdParams;
	int m_timeScale; //
	int m_steps; //# steps
  float m_time;
	int m_stepMax;
	bool m_exeTrajFlag; //true for executing
	void InitTrajGeneration(STrajInitParams *inputParams, int time);  	
	void ExeTrajGeneration(void);
  void FinishTrajGeneration(void);
	void RotationMatrix(void);
	double* CrossProduct(double* a, double* b);
  void RotationToRPYT(double x[], double y[], double z[]);
	//test function	
	void PrintLog(void);
     
	 
};

/*****************************
Constructor of CTrajGeneration
******************************/
CTrajGeneration::CTrajGeneration(ros::NodeHandle &n)
{
  pub = n.advertise<geometry_msgs::Twist>("/cmd_vel",10); 
  baseCmd.linear.x = baseCmd.linear.y = baseCmd.linear.z = 0;
  baseCmd.angular.x = baseCmd.angular.y = baseCmd.angular.z = 0;
  sFlightCmdParams.roll = sFlightCmdParams.pitch = sFlightCmdParams.yaw = sFlightCmdParams.throttle = 0;
  m_samplingTime = 0.01;//unit:sec
  m_steps = 0;
  m_timeScale = 0;
  m_time = 0.0;
  m_stepMax = 0;
  m_exeTrajFlag = false;
  
  for (int i = 0; i < 4; i++)
  {
    sTrajParams.qt[i] = 0;
    sTrajParams.qdt[i] = 0;
    sTrajParams.qddt[i] = 0;
  }
}


/************************************
  Init polynomials of trajectories
*************************************/
void CTrajGeneration::InitTrajGeneration(STrajInitParams *inputParams, int time)
{
  STrajInitParams input =  *inputParams;
  m_timeScale = time;
  m_time = 0/m_timeScale; //if initial time is not zero, need to divided by timescal
  m_steps = 0;
  m_stepMax = m_timeScale/m_samplingTime;
  
  //Calculate parameters for four polynomials (x, y, z, yaw)
  for(int i = 0; i < 4; i++)
  {
    /*
    % compute the polynomial coefficients
    A = 6*(q1 - q0) - 3*(qd1+qd0)*tscal;
    B = -15*(q1 - q0) + (8*qd0 + 7*qd1)*tscal;
    C = 10*(q1 - q0) - (6*qd0 + 4*qd1)*tscal;
    E = qd0*tscal; % as the t vector has been normalized
    F = q0;
    */ 
      
    sPolyParams.A[i] = 6*(input.goalPos[i] - input.initPos[i]) - 3*(input.goalVel[i] + input.initVel[i])*m_timeScale;
    sPolyParams.B[i] = -15*(input.goalPos[i] - input.initPos[i]) + (8*input.initVel[i] + 7*input.goalVel[i])*m_timeScale;
    sPolyParams.C[i] = 10*(input.goalPos[i] - input.initPos[i]) - (6*input.initVel[i] + 4*input.goalVel[i])*m_timeScale;
    sPolyParams.E[i] = input.initVel[i]*m_timeScale;
    sPolyParams.F[i] = input.initPos[i];
  }
  m_exeTrajFlag = true;
}

/****************************
  Execute
*****************************/
void CTrajGeneration::ExeTrajGeneration( void )
{
  double t5, t4, t3, t2, t;
  t = m_time;
  t5 = pow(t,5);
  t4 = pow(t,4);
  t3 = pow(t,3);
  t2 = pow(t,2);

  for(int i = 0; i < 4; i++)
  {
    sTrajParams.qt[i] = sPolyParams.A[i]*t5 + sPolyParams.B[i]*t4 + sPolyParams.C[i]*t3 + sPolyParams.E[i]*t + sPolyParams.F[i];     
    sTrajParams.qdt[i] = (5*sPolyParams.A[i]*t4 + 4*sPolyParams.B[i]*t3 + 3*sPolyParams.C[i]*t2 + sPolyParams.E[i])/m_timeScale;     
    sTrajParams.qddt[i] = ((20*sPolyParams.A[i]*t3 + 12*sPolyParams.B[i]*t2 + 6*sPolyParams.C[i]*t)/m_timeScale)/m_timeScale;     
  }
  //gLog[m_steps] = sTrajParams.qdt[0];
  //std::cout << "x" << " " << sTrajParams.qt[0] << std::endl;
  //std::cout << "y" << " " << sTrajParams.qt[1] << std::endl;
  //std::cout << "z" << " " << sTrajParams.qdt[2] << std::endl;
  //RotationMatrix();
  
  baseCmd.linear.x = sTrajParams.qdt[0];
  baseCmd.linear.y = sTrajParams.qdt[1];
  baseCmd.linear.z = sTrajParams.qdt[2];

  m_time = m_time + m_samplingTime/m_timeScale ;
  m_steps++;
}

/************
  Finish
************/
void CTrajGeneration::FinishTrajGeneration(void)
{
   m_steps = 0;
   m_time = 0;
   m_flightCmd.roll = 0;
   m_flightCmd.pitch = 0;
   m_flightCmd.yaw = 0;
   m_flightCmd.throttle = GRAVITY*MASS_OF_QUAD;
   //m_stepMax = 0;
}


/*****************************
  Calculate Rotation Matrix
******************************/
void CTrajGeneration::RotationMatrix( void )
{
  double t[3], x_c[3]= {0};
  double z_b[3];
  t[0] = sTrajParams.qddt[0];
  t[1] = sTrajParams.qddt[1];
  t[2] = sTrajParams.qddt[2] + GRAVITY;
  double tnorm = sqrt(pow(t[0],2) + pow(t[1],2) + pow(t[2],2));
  //caculate zb
  z_b[0] = t[0]/tnorm;
  z_b[1] = t[1]/tnorm;
  z_b[2] = t[2]/tnorm;
  
  x_c[0] = cos(sTrajParams.qt[3]);
  x_c[1] = sin(sTrajParams.qt[3]);
  x_c[2] = 0;

  //calculate yb
  double* temp;
  double y_b[3];
  temp = CrossProduct(z_b, x_c);
  double ynorm = sqrt(pow(temp[0],2) + pow(temp[1],2) + pow(temp[2],2));
  y_b[0] = temp[0]/ynorm;
  y_b[1] = temp[1]/ynorm;
  y_b[2] = temp[2]/ynorm;
  
  //caculate xb
  double x_b[3];
  temp = CrossProduct(y_b, z_b);
  x_b[0] = temp[0];
  x_b[1] = temp[1];
  x_b[2] = temp[2];
  

  RotationToRPYT(x_b,y_b,z_b);
  sFlightCmdParams.throttle = MASS_OF_QUAD*sqrt(pow(t[0],2) + pow(t[1],2) + pow(t[2],2));
  
  m_flightCmd.roll = sFlightCmdParams.roll;
  m_flightCmd.pitch = sFlightCmdParams.pitch;
  m_flightCmd.yaw = sFlightCmdParams.yaw;
  m_flightCmd.throttle = sFlightCmdParams.throttle;
  /* debug roll, pitch, yaw
  std::cout << "roll:" << sFlightCmdParams.roll <<\
  " " << "pitch:" << sFlightCmdParams.pitch << " "<<\ 
  "yaw:" << sFlightCmdParams.yaw << " "<< \
  "throttle:" << sFlightCmdParams.throttle << std::endl;
  */
  /* debug Rotation Matrix
  std::cout << x_b[0] << " " << y_b[0] << " "<<z_b[0] << std::endl;
  std::cout << x_b[1] << " " << y_b[1] << " "<<z_b[1] << std::endl;
  std::cout << x_b[2] << " " << y_b[2] << " "<<z_b[2] << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  */
  }

/*************
  Traj loop
**************/
void CTrajGeneration::TrajLoop( void )
{
   if ((m_steps < m_stepMax) && (m_exeTrajFlag == true))
   {
     ExeTrajGeneration();
     pub.publish(baseCmd);
   }
   
   if( m_steps == m_stepMax )
   {
     m_exeTrajFlag = false;
     FinishTrajGeneration();
     pub.publish(baseCmd);
     //PrintLog();
   }
}

/***************
  Cross Product
****************/
double* CTrajGeneration::CrossProduct(double* a, double* b)
{
  static double output[3] = {0};
  output[0] = a[1]*b[2] - a[2]*b[1];
  output[1] = a[2]*b[0] - a[0]*b[2];
  output[2] = a[0]*b[1] - a[1]*b[0];
  
  return output;
}

/*************************************
Convert Rotation Matrix to Flight Cmds
**************************************/
void CTrajGeneration::RotationToRPYT(double x[], double y[], double z[])
{
  double eps = 2.22e-15;
  double sp,cp;
  /*
  std::cout << "z[1] = "  << std::abs(z[1])<< std::endl;
  std::cout << "z[2] = "  << abs(z[2])<< std::endl;
  std::cout << "eps = "  << eps << std::endl;
  */
  if( (std::abs(z[2])<eps) && (std::abs(z[1]) < eps) )
  {
    std::cout << "what " << std::endl;
    sFlightCmdParams.roll = 0;
    sFlightCmdParams.pitch =  atan2(z[0],z[2]);
    sFlightCmdParams.yaw =  atan2(x[1],y[1]);
  }
  else
  {
    sFlightCmdParams.roll =  atan2(-z[1],z[2]);
    sp = sin(sFlightCmdParams.roll);
    cp = cos(sFlightCmdParams.roll);
    sFlightCmdParams.pitch = atan2( z[0], cp*z[2] - sp*z[1]);
    sFlightCmdParams.yaw = atan2( -y[0], x[0]);
  }
  sFlightCmdParams.roll = sFlightCmdParams.roll*180/PI; 
  sFlightCmdParams.pitch = sFlightCmdParams.pitch*180/PI; 
  sFlightCmdParams.yaw = sFlightCmdParams.yaw*180/PI; 
}


/***********************
  Print Log for Debug
************************/
void CTrajGeneration::PrintLog(void)
{
  for( int i = 0; i < 100; i++ )
  {
    std::cout << "qt[" << i <<"] = "  << gLog[i]<< std::endl;
  }
}


/*************
  TestFuntion
**************/
bool CTrajGeneration::TestFunction( void )
{
  //Cross Product Test
  /*
  double a[3] = {-0.0061, 0, 1};
  double b[3] = {1, 0, 0};
  double* c; 
  c = CrossProduct(a,b);
  std::cout << "c[0] = "  << c[0]<< std::endl;
  std::cout << "c[1] = "  << c[1]<< std::endl;
  std::cout << "c[2] = "  << c[2]<< std::endl;
  */

  
  //Temp entry point
  STrajInitParams inputParams;
  for (int i = 0; i < 4; i++)
  {
    inputParams.initPos[i] = 0;
    inputParams.initVel[i] = 0;
    inputParams.goalPos[i] = 0;
    inputParams.goalVel[i] = 0;
  }
  inputParams.initPos[0] = 0;
  inputParams.initPos[1] = 0;
  inputParams.initPos[2] = 0;
  inputParams.initPos[3] = 0;
  inputParams.goalPos[0] = 5;
  inputParams.goalPos[1] = 5;
  inputParams.goalPos[2] = 5;
  inputParams.goalPos[3] = 0;

  int time = 10;
  InitTrajGeneration(&inputParams, time);   
  

 
  //Gazebo Test 
  /*
  static int i = 0;
  i = i + 1;
  if( i>3)
  {
    i = 0;
  }
  //baseCmd.linear.z = baseCmd.linear.z + 0.01*i;
  //baseCmd.linear.z = 0;
  baseCmd.linear.x = x[i];
  pub.publish(baseCmd);
  */
}


/**********************************************************
  Add a Callback function Subscribe Prediction Algorithm
***********************************************************/

/******
  main
*******/
int main(int argc, char **argv)
{
  ros::init(argc, argv, "test1Fly");
  ros::NodeHandle nh;
  CTrajGeneration driver(nh);
  //ros::Rate loop_rate(0.5);
  ros::Rate loop_rate(100);
  
  driver.TestFunction();

  while(ros::ok())
  {
   //std::cout << "test" << std::endl;
   //driver.TestFunction();
   
   
   driver.TrajLoop();


   ros::spinOnce();
   loop_rate.sleep();
  }

  return 0;
}







