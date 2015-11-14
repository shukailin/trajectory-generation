cmake_minimum_required(VERSION 2.8.3)
project(traj_generation)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  geometry_msgs
  message_generation
)

add_message_files(
  FILES
  FlightCmd.msg
)

generate_messages(
  DEPENDENCIES
  std_msgs
)

catkin_package(CATKIN_DEPENDS roscpp rospy std_msgs message_runtime)

###########
## Build ##
###########

include_directories(${catkin_INCLUDE_DIRS})
add_executable(test1 src/test1.cpp)
target_link_libraries(test1 ${catkin_LIBRARIES})
add_dependencies(test1 traj_generation_generate_message_cpp)

add_executable(test1Fly src/test1Fly.cpp)
target_link_libraries(test1Fly ${catkin_LIBRARIES})
add_dependencies(test1Fly traj_generation_generate_message_cpp)

add_executable(trajGeneration src/trajGeneration.cpp)
target_link_libraries(trajGeneration ${catkin_LIBRARIES})
add_dependencies(trajGeneration traj_generation_generate_message_cpp)

add_executable(flightControl src/flightControl.cpp)
target_link_libraries(flightControl ${catkin_LIBRARIES})
add_dependencies(flightControl traj_generation_generate_message_cpp)


