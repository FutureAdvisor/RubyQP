require 'bundler/gem_tasks'

#
# Because this is a private gem, remove the release task so it doesn't
# accidentally get released by accident.
#

Rake::TaskManager.class_eval do
  def remove_task(task_name)
    @tasks.delete(task_name.to_s)
  end
end

Rake.application.remove_task 'release'
