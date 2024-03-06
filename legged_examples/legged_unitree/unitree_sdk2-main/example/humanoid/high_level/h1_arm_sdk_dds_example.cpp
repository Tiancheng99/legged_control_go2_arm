#include <array>
#include <chrono>
#include <iostream>
#include <thread>

#include <unitree/idl/go2/LowCmd_.hpp>
#include <unitree/robot/channel/channel_publisher.hpp>

static const std::string kTopicArmSDK = "rt/arm_sdk";
constexpr float kPi = 3.141592654;
constexpr float kPi_2 = 1.57079632;

enum JointIndex {
  // Right leg
  kRightHipYaw = 8,
  kRightHipRoll = 0,
  kRightHipPitch = 1,
  kRightKnee = 2,
  kRightAnkle = 11,
  // Left leg
  kLeftHipYaw = 7,
  kLeftHipRoll = 3,
  kLeftHipPitch = 4,
  kLeftKnee = 5,
  kLeftAnkle = 10,

  kWaistYaw = 6,

  kNotUsedJoint = 9,

  // Right arm
  kRightShoulderPitch = 12,
  kRightShoulderRoll = 13,
  kRightShoulderYaw = 14,
  kRightElbow = 15,
  // Left arm
  kLeftShoulderPitch = 16,
  kLeftShoulderRoll = 17,
  kLeftShoulderYaw = 18,
  kLeftElbow = 19,

};

int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " networkInterface" << std::endl;
    exit(-1);
  }

  unitree::robot::ChannelFactory::Instance()->Init(0, argv[1]);

  unitree::robot::ChannelPublisherPtr<unitree_go::msg::dds_::LowCmd_>
      arm_sdk_publisher;
  unitree_go::msg::dds_::LowCmd_ msg;

  arm_sdk_publisher.reset(
      new unitree::robot::ChannelPublisher<unitree_go::msg::dds_::LowCmd_>(
          kTopicArmSDK));
  arm_sdk_publisher->InitChannel();

  std::array<JointIndex, 8> arm_joints = {
      JointIndex::kLeftShoulderPitch,  JointIndex::kLeftShoulderRoll,
      JointIndex::kLeftShoulderYaw,    JointIndex::kLeftElbow,
      JointIndex::kRightShoulderPitch, JointIndex::kRightShoulderRoll,
      JointIndex::kRightShoulderYaw,   JointIndex::kRightElbow};

  float weight = 0.f;
  float weight_rate = 0.2f;

  float kp = 60.f;
  float kd = 1.5f;
  float dq = 0.f;
  float tau_ff = 0.f;

  float control_dt = 0.02f;
  float max_joint_velocity = 0.5f;

  float delta_weight = weight_rate * control_dt;
  float max_joint_delta = max_joint_velocity * control_dt;
  auto sleep_time =
      std::chrono::milliseconds(static_cast<int>(control_dt / 0.001f));

  std::array<float, 8> init_pos{};

  std::array<float, 8> target_pos = {0.f, kPi_2,  0.f, kPi_2,
                                     0.f, -kPi_2, 0.f, kPi_2};

  // wait for init
  std::cout << "Press ENTER to init arms ...";
  std::cin.get();

  // set init pos
  std::cout << "Initailizing arms ...";
  float init_time = 5.0f;
  int init_time_steps = static_cast<int>(init_time / control_dt);

  for (int i = 0; i < init_time_steps; ++i) {
    // increase weight
    weight += delta_weight;
    weight = std::clamp(weight, 0.f, 1.f);
    std::cout << weight << std::endl;

    // set weight
    msg.motor_cmd().at(JointIndex::kNotUsedJoint).q(weight * weight);

    // set control joints
    for (int j = 0; j < init_pos.size(); ++j) {
      msg.motor_cmd().at(arm_joints.at(j)).q(init_pos.at(j));
      msg.motor_cmd().at(arm_joints.at(j)).dq(dq);
      msg.motor_cmd().at(arm_joints.at(j)).kp(kp);
      msg.motor_cmd().at(arm_joints.at(j)).kd(kd);
      msg.motor_cmd().at(arm_joints.at(j)).tau(tau_ff);
    }

    // send dds msg
    arm_sdk_publisher->Write(msg);

    // sleep
    std::this_thread::sleep_for(sleep_time);
  }

  std::cout << "Done!" << std::endl;

  // wait for control
  std::cout << "Press ENTER to start arm ctrl ..." << std::endl;
  std::cin.get();

  // start control
  std::cout << "Start arm ctrl!" << std::endl;
  float period = 5.f;
  int num_time_steps = static_cast<int>(period / control_dt);

  std::array<float, 8> current_jpos_des{};

  // lift arms up
  for (int i = 0; i < num_time_steps; ++i) {
    // update jpos des
    for (int j = 0; j < init_pos.size(); ++j) {
      current_jpos_des.at(j) +=
          std::clamp(target_pos.at(j) - current_jpos_des.at(j),
                     -max_joint_delta, max_joint_delta);
    }

    // set control joints
    for (int j = 0; j < init_pos.size(); ++j) {
      msg.motor_cmd().at(arm_joints.at(j)).q(current_jpos_des.at(j));
      msg.motor_cmd().at(arm_joints.at(j)).dq(dq);
      msg.motor_cmd().at(arm_joints.at(j)).kp(kp);
      msg.motor_cmd().at(arm_joints.at(j)).kd(kd);
      msg.motor_cmd().at(arm_joints.at(j)).tau(tau_ff);
    }

    // send dds msg
    arm_sdk_publisher->Write(msg);

    // sleep
    std::this_thread::sleep_for(sleep_time);
  }

  // put arms down
  for (int i = 0; i < num_time_steps; ++i) {
    // update jpos des
    for (int j = 0; j < init_pos.size(); ++j) {
      current_jpos_des.at(j) +=
          std::clamp(init_pos.at(j) - current_jpos_des.at(j), -max_joint_delta,
                     max_joint_delta);
    }

    // set control joints
    for (int j = 0; j < init_pos.size(); ++j) {
      msg.motor_cmd().at(arm_joints.at(j)).q(current_jpos_des.at(j));
      msg.motor_cmd().at(arm_joints.at(j)).dq(dq);
      msg.motor_cmd().at(arm_joints.at(j)).kp(kp);
      msg.motor_cmd().at(arm_joints.at(j)).kd(kd);
      msg.motor_cmd().at(arm_joints.at(j)).tau(tau_ff);
    }

    // send dds msg
    arm_sdk_publisher->Write(msg);

    // sleep
    std::this_thread::sleep_for(sleep_time);
  }

  // stop control
  std::cout << "Stoping arm ctrl ...";
  float stop_time = 2.0f;
  int stop_time_steps = static_cast<int>(stop_time / control_dt);

  for (int i = 0; i < stop_time_steps; ++i) {
    // increase weight
    weight -= delta_weight;
    weight = std::clamp(weight, 0.f, 1.f);

    // set weight
    msg.motor_cmd().at(JointIndex::kNotUsedJoint).q(weight);

    // send dds msg
    arm_sdk_publisher->Write(msg);

    // sleep
    std::this_thread::sleep_for(sleep_time);
  }

  std::cout << "Done!" << std::endl;

  return 0;
}
