#include <iostream>
#include <ros/ros.h>
#include <geometry_msgs/Twist.h>

class CRobotDriver
{
  public:
	CRobotDriver(ros::NodeHandle &nh)
	{
	  n = nh;
	  pub = n.advertise<geometry_msgs::Twist>("/cmd_vel_mux/input/teleop",10); 
	  baseCmd.linear.x = baseCmd.linear.y = baseCmd.linear.z = 0;
	  baseCmd.angular.x = baseCmd.angular.y = baseCmd.angular.z = 0;
	}
	bool TestFunction(void);	  	

  private:
	ros::NodeHandle n;
	ros::Publisher pub;
	geometry_msgs::Twist baseCmd;
};

bool CRobotDriver::TestFunction(void)
{
  static float i = 0;
  i = i + 1;
  //baseCmd.linear.z = baseCmd.linear.z + 0.01*i;
  baseCmd.linear.x = 0.5;
  pub.publish(baseCmd);
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "test1");
  ros::NodeHandle nh;
  CRobotDriver driver(nh);
  ros::Rate loop_rate(1000);
  while(ros::ok())
  {
   driver.TestFunction();
   ros::spinOnce();
   loop_rate.sleep();
  }

  return 0;
}


