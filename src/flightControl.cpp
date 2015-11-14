#include <iostream>
#include <ros/ros.h>
#include <traj_generation/FlightCmd.h>

class CFlightControl
{
  public:
  	CFlightControl(ros::NodeHandle &nh);     
  private:
	ros::Subscriber sub;
	void FlightCallback(const traj_generation::FlightCmd& msg);
};


void CFlightControl::FlightCallback(const traj_generation::FlightCmd& msg)
{

}


CFlightControl::CFlightControl(ros::NodeHandle &n)
{
   sub = n.subscribe("flightCmd", 10, &CFlightControl::FlightCallback, this);
}


int main(int argc, char **argv)
{
  ros::init(argc, argv, "flightControl");
  ros::NodeHandle nh;
  CFlightControl CflightControl(nh);
  
  ros::Rate loop_rate(100);
  
  while(ros::ok())
  {
   ros::spinOnce();
   loop_rate.sleep();
  }
  return 0;
}


